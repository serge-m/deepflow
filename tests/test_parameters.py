#!/usr/bin/python

import unittest
from py_fastdeepflow import fastdeepflow
import cv2
import os
import numpy
import ctypes

path_dir_data = os.path.join(os.path.dirname(__file__), 'data/')


def imread(path):
    """Read image from file

    param: path path to source image

    return: image as numpy array
    """
    img = cv2.imread(path)

    if img is None:
        raise Exception("Unable to open '{}'".format(path))

    if img.ndim != 1 and img.ndim != 3:
        raise Exception("Unsupported number of dimensions")

    if img.ndim == 3 and img.shape[-1] == 3:
        img = img[:,:,::-1]

    return img


class TestParameters(unittest.TestCase):
    def prepare(self, dir_name, path_flow_reference):

        img0 = imread(os.path.join(path_dir_data, 'artificial/%s/img0.png' % dir_name))
        img1 = imread(os.path.join(path_dir_data, 'artificial/%s/img1.png' % dir_name))

        flow_reference = fastdeepflow.read_flow(os.path.join(path_dir_data, ('artificial/%s/' + path_flow_reference) % dir_name))

        return img0, img1, flow_reference

    def test_middleburry(self):
        img0, img1, flow_reference = self.prepare(dir_name='48_0', path_flow_reference='48_midl.flo')

        params = fastdeepflow.optical_flow_params_t()
        fastdeepflow.lib.optical_flow_params_middlebury(ctypes.byref(params))
        flow = fastdeepflow.calc_flow(img0, img1, params)

        self.assertTrue(numpy.array_equal(numpy.array(flow), numpy.array(flow_reference)))

    def test_middleburry__min_size_10(self):
        img0, img1, flow_reference = self.prepare(dir_name='48_0', path_flow_reference='48_midl_minsize10.flo')

        params = fastdeepflow.optical_flow_params_t()
        fastdeepflow.lib.optical_flow_params_middlebury(ctypes.byref(params))
        params.min_size = 10
        flow = fastdeepflow.calc_flow(img0, img1, params)

        self.assertTrue(numpy.array_equal(numpy.array(flow), numpy.array(flow_reference)))


if __name__ == '__main__':
    unittest.main()
