#!/usr/bin/env python
from __future__ import print_function
import json
import logging
import os
import subprocess
import unittest

import dxpy
import yaml


TEST_DIRECTORY = os.path.dirname(os.path.realpath(__file__))
TEST_DATA_YAML = os.path.join(TEST_DIRECTORY, '../../../test_data.yml')
MAIN_APP_DIR = os.path.join(TEST_DIRECTORY, '../')
BUILD_CMD = ['dx', 'build', '-f', MAIN_APP_DIR]
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


def build_applet():
    applet_config = subprocess.check_output(BUILD_CMD).strip()
    return json.loads(applet_config)['id']


def remove_applet(applet_id):
    subprocess.check_call(['dx', 'rm', applet_id])


def remove_output(output_value):
    if isinstance(output_value, list):
        map(remove_output, output_value)
    elif isinstance(output_value, dict):
        map(remove_output, output_value.values())
    elif dxpy.is_dxlink(output_value):
        dxpy.remove(output_value)
    elif True: # TODO: test whether output_value is a job
        remove_output(output_value.describe()['output'])
    else:
        LOG.warn('Unrecognized output type: %s', output_value)


class VcfAnnotateTests(unittest.TestCase):
    def setUp(self):
        LOG.info('Loading test data')
        self.test_data = load_test_data()

        LOG.info('Building applets')
        self.applet_id = build_applet()

        self.jobs = []

    def tearDown(self):
        LOG.info('Removing applet')
        remove_applet(self.applet_id)

        LOG.info('Removing job outputs')
        remove_output(self.jobs)
        self.jobs = None

    def test_applet(self):
        self.assertIsNotNone(self.test_data)
        self.assertIsNotNone(self.applet_id)

        applet = dxpy.DXApplet(self.applet_id)

        job_input = {
            'variants_vcfgzs': [
                dxpy.dxlink(self.test_data['hg001_wes_vcfgz']),
                dxpy.dxlink(self.test_data['hg002_wes_vcfgz'])
            ],
            'variants_vcfgztbis': [
                dxpy.dxlink(self.test_data['hg001_wes_vcfgz_tbi']),
                dxpy.dxlink(self.test_data['hg002_wes_vcfgz_tbi'])
            ],
            'build_version': 'b38'
        }
        self.jobs.append(applet.run(job_input))
        LOG.info('Running full test in job: %d.', self.jobs[-1].get_id())

        for job in self.jobs:
            # Raises an error if the job fails or does not complete before the timeout
            job.wait_on_done()


if __name__ == '__main__':
    unittest.main()
