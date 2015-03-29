#!/usr/bin/python

import unittest
from py_fastdeepflow import fastdeepflow
import cv2
import os
import numpy
import ctypes

path_dir_data = os.path.join(os.path.dirname(__file__), 'data/')

class TestStringMethods(unittest.TestCase):


    def test_image_reading(self):
        img = fastdeepflow.lib.color_image_load(path_dir_data + 'sintel1.png')
        i = fastdeepflow.numpy_from_color_image_t(img.contents)
        self.assertTrue(i.shape == (436, 1024, 3))

    def test1(self):
        img0 = cv2.imread(path_dir_data + 'sintel1.png')[:,:,::-1]
        img1 = cv2.imread(path_dir_data + 'sintel2.png')[:,:,::-1]

    def test_params(self):
        params = fastdeepflow.optical_flow_params_t()
        fastdeepflow.lib.optical_flow_params_default(ctypes.byref(params))
        self.assertTrue(params.n_solver_iteration==25)
        self.assertTrue(abs(params.sor_omega-1.6) < 1e-5)

    def test_of(self):
        img0 = cv2.imread(path_dir_data + 'sintel1.png')[:,:,::-1]
        img1 = cv2.imread(path_dir_data + 'sintel2.png')[:,:,::-1]
        u, v = fastdeepflow.calc_flow(img0, img1)
        u1, v1 = fastdeepflow.read_flow(path_dir_data + 'sintel_cmd.flo')
        self.assertTrue(numpy.allclose(u, u1, atol=1e-8))
        self.assertTrue(numpy.allclose(v, v1, atol=1e-8))

    def test_wrap(self):
        img0 = cv2.imread(path_dir_data + 'sintel1.png')[:,:,::-1]
        img1 = cv2.imread(path_dir_data + 'sintel2.png')[:,:,::-1]
        u, v = fastdeepflow.calc_flow(img0, img1)
        img0warped = fastdeepflow.warp_image(img0, u, v)
        img1warped = fastdeepflow.warp_image(img1, u, v)
        
        d01 = numpy.sum(numpy.abs(img0-img1.astype('float32')))
        d0w1= numpy.sum(numpy.abs(img0warped-img1))
        d01w= numpy.sum(numpy.abs(img0-img1warped))
        print "d01 {}, d0w1 {}, d01w {}".format(d01, d0w1, d01w)
        self.assertTrue(d01w<d01)
        self.assertTrue(d0w1>d01)
    

if __name__ == '__main__':
    unittest.main()
