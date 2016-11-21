#!/bin/bash

function usage(){
echo \
"Usage: $(basename $0) [-o org] [-p project] [OPTIONS]

Displays the current permissions for members in an organization, as well
as the permissions of non-members in projects billed to this organization.

Options:
	-o Organization to show permissions for
	-p Restrict output to these projects (may be given more than once)
	-m Show member permissions only; no projects
	-e Show permissions in projects eleveated above a user's project-access
	-s Show permissions separately for each project
    
NOTE: You must have already logged into DNAnexus and have either ADMIN permission
on the org.
"

exit 1
}


ORG=""
PROJECTS=()
MEMBER_ONLY=0
ELEVATED=0
SEPARATE_PERMS=0

while getopts o:p:mesh flag; do
	case $flag in
	o)
		ORG="$OPTARG";
		;;
	p)
		PROJECTS+=("$OPTARG");
		;;
	m)
		MEMBER_ONLY=1
		;;
	e)
		ELEVATED=1
		;; 
	s)
	    SEPARATE_PERMS=1
	    ;;
	h)
		usage
		;;
	\?)
		echo "Invalid option: -$OPTARG"
		usage
		;;
	esac
done

PROJ_ADMINS_JS=$(mktemp)
PROJ_MEMBERS_JS=$(mktemp)

# first, find all of the admins
echo "Adminstrators of $ORG:"
echo -e "Name\tID"
dx find org members $ORG --level ADMIN --json | tee $PROJ_ADMINS_JS | jq -r '.[] | [.describe.first, .describe.last, .id] | @csv' | sed 's/,/ /' | tr ',' '\t' | sed 's/"//g'

echo ""
echo ""
echo "Members of $ORG:"
echo -e "Name\tID\tAccess\tBillable?"
dx find org members $ORG --level MEMBER --json | tee $PROJ_MEMBERS_JS | jq -r '.[] | [.describe.first, .describe.last, .id, .projectAccess, .allowBillableActivities] | @csv' | sed 's/,/ /' | tr ',' '\t' | sed 's/"//g'

echo ""
echo ""

# convert levels to int, lower int is "higher" permissions
declare -A PERM_LEVEL
PERM_LEVEL["ADMINISTER"]=1
PERM_LEVEL["CONTRIBUTE"]=2
PERM_LEVEL["UPLOAD"]=3
PERM_LEVEL["VIEW"]=4

# Get an array of users -> project access level in the org
declare -A ORG_LEVEL

# first, admins always have "ADMINSTER" in the org
while read n; do
    ORG_LEVEL["$n"]="ADMINISTER"
done < <(cat $PROJ_ADMINS_JS | jq -r '.[] | [.id] | @csv' | sed 's/"//g')

# Now, go through the members and give them the actual access
while read n; do
    NAME=$(echo "$n" | cut -d, -f1)
    LEVEL=$(echo "$n" |cut -d, -f2)
    ORG_LEVEL["$NAME"]="$LEVEL"
done < <(cat $PROJ_MEMBERS_JS | jq -r '.[] | [.id, .projectAccess] | @csv' | sed 's/"//g')

rm $PROJ_ADMINS_JS
rm $PROJ_MEMBERS_JS

# Get a list of projects
PROJ_LIST_FN=$(mktemp)

#set -x

if test "$PROJECTS"; then
    for p in "${!PROJECTS[@]}"; do
        echo "${PROJECTS[$p]}" >> $PROJ_LIST_FN
    done
else
	# go through the org to get the list of projects
	PROJ_START_STR=" "
	while test "$PROJ_START_STR"; do
		PROJ_IN_JSON=$(echo "'"'{"limit": 1000, "describe" : true'$PROJ_START_STR"}'")
		PROJ_RES_JSON=$(eval dx api $ORG findProjects "$PROJ_IN_JSON")
		
		echo $PROJ_RES_JSON | jq -r '.results | .[] | [.id] | @csv' | sed 's/"//g' >> $PROJ_LIST_FN
		
		NEXT_PROJ=$(echo $PROJ_RES_JSON | jq .next)
		if test ! "$NEXT_PROJ" == "null"; then
			PROJ_START_STR=", \"starting\" : $NEXT_PROJ"
		else
			PROJ_START_STR=""
		fi
	done
fi

SEP_VAL="\t"
HDR_STR="Projects in $ORG:\nName\tID\tCreator"
if test $SEPARATE_PERMS -ne 0; then
    SEP_VAL=": "
    HDR_STR="Project Permissions in $ORG:"
fi

ALL_PERMS_FN=$(mktemp)

echo -e "$HDR_STR"
while read p; do
    PROJ_JSON=$(dx describe --json --verbose $p)
    PROJ_NM=$(echo "$PROJ_JSON" | jq -r .name)
    PROJ_ID=$(echo "$PROJ_JSON" | jq -r .id)
    PROJ_CREATOR=$(echo "$PROJ_JSON" | jq -r '.createdBy["user"]')
    PROJ_NAME="${PROJ_NM} ($PROJ_ID)"
    if test $SEPARATE_PERMS -eq 0; then
        PROJ_NAME="${PROJ_NM}\t${PROJ_ID}"
    fi
    
    echo -e "${PROJ_NAME}${SEP_VAL}${PROJ_CREATOR}"

    # get a list of user -> permission level in the project
    unset PROJ_LEVEL
    declare -A PROJ_LEVEL
    
    for u in $(echo "$PROJ_JSON" | jq -r '.permissions|keys[]'); do
#        echo "U = $u"
        PROJ_LEVEL["$u"]=$(echo "$PROJ_JSON" | jq -r ".permissions[\"$u\"]")
    done
    
    unset PROJ_PRINT
    declare -A PROJ_PRINT
    # Now, iterate over proj_level and only keep the "relevant" ones
    for k in "${!PROJ_LEVEL[@]}"; do
#        echo "proj_level $k ${PROJ_LEVEL[$k]}"
    
        # get the level in the org
        USER_ORG_LEVEL="${ORG_LEVEL[$k]}"
        
        # if not present in the org, OR if (ELEVATED is True AND USERORG_LEVEL is higher than the given level), add to PROJ_PRINT
        if [[ -z "$USER_ORG_LEVEL" ]] || ( [[ $ELEVATED -ne 0 ]] && [[ ${PERM_LEVEL[$USER_ORG_LEVEL]} -gt ${PERM_LEVEL[${PROJ_LEVEL[$k]}]} ]] ) ; then
            PROJ_PRINT["$k"]=${PROJ_LEVEL[$k]}
        fi
    done
    
    PROJ_PERMS_FN=$(mktemp)
    
    # OK, now print the user levels
    for k in "${!PROJ_PRINT[@]}"; do
        if test $SEPARATE_PERMS -eq 0; then
            echo -en "$PROJ_NAME\t"
        fi
        echo -e "$k\t${PROJ_PRINT[$k]}"
    done > $PROJ_PERMS_FN
    
    if test $SEPARATE_PERMS -ne 0; then
        cat $PROJ_PERMS_FN
        echo ""
    else
        cat $PROJ_PERMS_FN >> $ALL_PERMS_FN
    fi  

done < $PROJ_LIST_FN

rm $PROJ_LIST_FN

echo ""

if test $SEPARATE_PERMS -eq 0; then
    echo "Project Permissions in $ORG:"
    echo -e "Project Name\tID\tUser\tPermission Level"
    cat $ALL_PERMS_FN
    echo ""
fi
