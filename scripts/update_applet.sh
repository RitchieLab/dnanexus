#!/bin/bash

#set -x

# First find the applet we're replacing
CURR_APPLET=$(dx find data --class applet --name "$1" --brief)
N_CURR_APPLETS=$(dx find data --class applet --name "$1" --brief | wc -l)

CURR_PROJECT=$(dx env --bash | grep 'DX_PROJECT_CONTEXT_ID' | sed 's/[^=]*=//')
CURR_DXDIR=$(dx env | grep 'Current folder' | sed 's/^Current folder[ \t]*//')

if test "$N_CURR_APPLETS" -gt 1; then
	echo "WARNING: Ambiguous applet, checking for current directory"
	CURR_APPLET=$(dx find data --class applet --name "$1" --project "$CURR_PROJECT" --folder "$CURR_DXDIR" --norecurse --brief)
fi
	
if test -z "$CURR_APPLET"; then
	echo "ERROR: Applet to update not found!"
	exit 1
fi

PROJECT_ID=$(echo $CURR_APPLET | sed 's/:.*//')
OLD_APPLET_ID=$(echo $CURR_APPLET | sed 's/[^:]*://')
APPLET_FOLDER=$(dx describe $CURR_APPLET | grep '^Folder' | sed 's/^Folder *//')

# Get a list of workflows containing the old applet
WORKFLOW_UPDATE_FN=$(mktemp)
for f in $(dx find data --class workflow --brief); do 
	if test $(dx describe $f | grep $OLD_APPLET_ID | wc -l) -gt 0; then 
		echo $f; 
	fi; 
done > $WORKFLOW_UPDATE_FN

ID_FN=$(mktemp)

dx build -a "$1" -d "$APPLET_FOLDER/" > $ID_FN

if test "$?" -ne 0; then
	echo "Error: Build unsuccessful, aborting!"
	exit 2
fi

# build the new applet, archiving the old one and putting the new one in the current location
NEW_APPLET_ID=$(cat $ID_FN | sed -e 's/.*: *"//' -e 's/"}$//')
rm $ID_FN

if test -z "$NEW_APPLET_ID"; then
	echo "ERROR: No new applet ID found, aborting!"
	exit 3
fi

# Now for every workflow found above, find and update the stage
for w in $(cat $WORKFLOW_UPDATE_FN); do
	echo "updating Workflow $(dx describe --name $w)"
	dx describe $w | grep -B 1 $OLD_APPLET_ID | head -1 | sed 's/.* \(stage-.*$\)/\1/' | xargs -n1 -i dx update stage $w {} --executable $NEW_APPLET_ID
done

# Now, let's delete the old applet so we can't use it!
dx rm $CURR_APPLET

rm $WORKFLOW_UPDATE_FN
	
