#!/usr/bin/env python
from __future__ import print_function
import json
import logging
import os
import subprocess
import unittest
import uuid

import dxpy
import yaml


APP_NAME = 'vcf_annotate'
TEST_DIRECTORY = os.path.dirname(os.path.realpath(__file__))
TEST_DATA_YAML = os.path.join(TEST_DIRECTORY, '../../../test_data.yml')
MAIN_APP_DIR = os.path.join(TEST_DIRECTORY, '../')
MKDIR_CMD = ['dx', 'mkdir', '-p']
BUILD_CMD = ['dx', 'build', MAIN_APP_DIR, '-a', '--destination']
LOG = logging.getLogger()


def load_test_data():
    with open(TEST_DATA_YAML, 'rU') as fh:
        test_data = yaml.safe_load(fh)
    if dxpy.APISERVER_HOST.startswith('staging') and 'staging' in test_data:
        host_data = test_data['staging']
    elif 'production' in test_data:
        host_data = test_data['production']
    else:
        raise AssertionError(
            'Test data YAML file does not contain a configuration for current host %s',
            dxpy.APISERVER_HOST
        )
    region = dxpy.describe(dxpy.PROJECT_CONTEXT_ID)['region']
    if region not in host_data:
        raise AssertionError(
            'Test data YAML file does not contain a configuration for the current '
            'region %s on the current host %s', dxpy.APISERVER_HOST, region
        )
    return host_data[region]


def build_applet(target_project):
    suffix = str(uuid.uuid4())
    unq_app_name = '{}_{}'.format(APP_NAME, suffix)
    target_relpath = '/{}/{}'.format(APP_NAME, unq_app_name)
    target_path = '{}:{}'.format(target_project, target_relpath)
    mkdir_cmd = MKDIR_CMD + [target_path]
    subprocess.check_output(mkdir_cmd)
    build_cmd = BUILD_CMD + ['{}/{}'.format(target_path, unq_app_name)]
    applet_config = subprocess.check_output(build_cmd).strip()
    applet_id = json.loads(applet_config)['id']
    return target_relpath, applet_id


def remove_applet(target_project, target_relpath):
    subprocess.check_call([
        'dx', 'rm', '-r', '{}:{}'.format(target_project, target_relpath)
    ])


def remove_output(output_value):
    if isinstance(output_value, list):
        map(remove_output, output_value)
    elif isinstance(output_value, dict):
        map(remove_output, output_value.values())
    elif dxpy.is_dxlink(output_value):
        dxpy.remove(output_value)
    elif True:  # TODO: test whether output_value is a job
        output = output_value.describe()['output']
        if output:
            remove_output(output)
    else:
        LOG.warn('Unrecognized output type: %s', output_value)


class VcfAnnotateTests(unittest.TestCase):
    @property
    def target_project(self):
        return self.test_data['project']

    def setUp(self):
        LOG.info('Loading test data')
        self.test_data = load_test_data()

        LOG.info('Building applets')
        self.target_relpath, self.applet_id = build_applet(self.target_project)

        self.jobs = []

    def tearDown(self):
        LOG.info('Removing applet')
        remove_applet(self.target_project, self.target_relpath)

        # TODO: this is probably no longer necessary
        LOG.info('Removing job outputs')
        remove_output(self.jobs)
        self.jobs = None

    def test_applet(self):
        self.assertIsNotNone(self.test_data)
        self.assertIsNotNone(self.applet_id)

        applet = dxpy.DXApplet(self.applet_id)

        job_input = {
            'variants_vcfgz': [
                dxpy.dxlink(self.test_data['hg001_wes_vcfgz']),
                dxpy.dxlink(self.test_data['hg002_wes_vcfgz'])
            ],
            'variants_vcfgztbi': [
                dxpy.dxlink(self.test_data['hg001_wes_vcfgz_tbi']),
                dxpy.dxlink(self.test_data['hg002_wes_vcfgz_tbi'])
            ],
            'clinvar_vcfgz': dxpy.dxlink(self.test_data['clinvar_vcfgz']),
            'clinvar_vcfgztbi': dxpy.dxlink(self.test_data['clinvar_vcfgztbi'])
        }
        self.jobs.append(applet.run(
            job_input, project=self.target_project, folder=self.target_relpath
        ))
        LOG.info('Running full test in job: %d.', self.jobs[-1].get_id())

        for job in self.jobs:
            # Raises an error if the job fails or does not complete before the timeout
            job.wait_on_done()


if __name__ == '__main__':
    unittest.main()
