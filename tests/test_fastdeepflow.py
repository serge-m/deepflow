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
    
    
    def test_wrap_grayscale(self):
        img0 = cv2.imread(path_dir_data + 'sintel1.png')[:,:,::-1]
        img1 = cv2.imread(path_dir_data + 'sintel2.png')[:,:,::-1]
        u, v = fastdeepflow.read_flow(path_dir_data + 'sintel_cmd.flo')
        img1warped = fastdeepflow.warp_image(img1, u, v)

        z = [fastdeepflow.warp_image(img1[:,:,k], u, v) for k in range(3)]
        zz = numpy.dstack(z, )
        
        
        self.assertTrue(numpy.allclose(img1warped, zz, atol=1e-8))


class ReadyFlow(unittest.TestCase):
    def setUp(self):
        self.u, self.v = fastdeepflow.read_flow(path_dir_data + 'sintel_cmd.flo')
        self.assertTrue(self.u.max()-self.u.min()>5)


class TestIO(ReadyFlow):
    def test_read_write(self):
        path_tmp_flow = os.path.join(path_dir_data, 'flow_read_write_test.flo')
        if os.path.exists(path_tmp_flow):
            os.remove(path_tmp_flow)
        self.assertFalse(os.path.exists(path_tmp_flow))

        self.assertEqual(self.u.shape, self.v.shape)

        fastdeepflow.write_flow(path_tmp_flow, self.u, self.v)
        self.assertTrue(os.path.exists(path_tmp_flow))
        u_loaded, v_loaded = fastdeepflow.read_flow(path_tmp_flow)
        self.assertTrue(numpy.array_equal(self.u, u_loaded))
        self.assertTrue(numpy.array_equal(self.v, v_loaded))


class TestExceptions(ReadyFlow):
    def test_read_write(self):
        with self.assertRaises(Exception):
            fastdeepflow.write_flow('123123/123/123123', self.u, self.v)

        with self.assertRaises(Exception):
            res = fastdeepflow.read_flow('123123/234234')


if __name__ == '__main__':
    unittest.main()
