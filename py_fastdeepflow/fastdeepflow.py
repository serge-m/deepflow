#!/usr/bin/python


import sys 
import os

import ctypes
import numpy

path_lib = os.path.join(os.path.dirname(__file__), 'libdeepflow.so')
lib = ctypes.cdll.LoadLibrary(path_lib)


class image_t(ctypes.Structure):
     _fields_ = [("width", ctypes.c_int),
                 ("height", ctypes.c_int),
                 ("stride", ctypes.c_int),
                 ("data", ctypes.POINTER(ctypes.c_float)),
                 ]

class color_image_t(ctypes.Structure):
     _fields_ = [("width", ctypes.c_int),
                 ("height", ctypes.c_int),
                 ("stride", ctypes.c_int),
                 ("c1", ctypes.POINTER(ctypes.c_float)),
                 ("c2", ctypes.POINTER(ctypes.c_float)),
                 ("c3", ctypes.POINTER(ctypes.c_float)),
                 ]


class optical_flow_params_t(ctypes.Structure):
    _fields_ = [("alpha", ctypes.c_float),
                ("beta", ctypes.c_float),
                ("gamma", ctypes.c_float),
                ("delta", ctypes.c_float),
                ("sigma", ctypes.c_float),
                ("bk", ctypes.c_float),
                ("eta", ctypes.c_float),
                ("min_size", ctypes.c_int),
                ("n_inner_iteration", ctypes.c_int),
                ("n_solver_iteration", ctypes.c_int),
                ("sor_omega", ctypes.c_float),
                ]

lib.readFlowFile.restype = ctypes.POINTER(ctypes.POINTER(image_t))
lib.readFlowFile.argtypes = [ctypes.c_char_p]
lib.color_image_load.restype = ctypes.POINTER(color_image_t)
lib.color_image_load.argtypes = [ctypes.c_char_p]
lib.image_delete.restype = None

lib.color_image_new.restype = ctypes.POINTER(color_image_t)
lib.image_new.restype = ctypes.POINTER(image_t)



def numpy_from_image_t(img_t):
    full = numpy.ctypeslib.as_array(img_t.data, shape=(img_t.height, img_t.stride))
    return full[:,:img_t.width].copy()

def numpy_from_color_image_t_1(img_t):
    shape=(img_t.height, img_t.stride)

    channels = [numpy.ctypeslib.as_array(c, shape=shape)[:,:img_t.width] for c in [img_t.c1, img_t.c2, img_t.c3]]

    return numpy.dstack(channels)


def numpy_from_color_image_t(img_t):
    shape=(3, img_t.height, img_t.stride)
    return numpy.ctypeslib.as_array(img_t.c1, shape=shape).transpose(1,2,0)[:,:,:img_t.width].copy()



def read_flow(path):
    t = lib.readFlowFile(path)

    u = numpy_from_image_t(t[0].contents)
    v = numpy_from_image_t(t[1].contents)

    lib.image_delete(t[0])
    lib.image_delete(t[1])

    return u, v


def fill_colot_image_t(img, new_c_img):
    h, w, s = new_c_img.height, new_c_img.width, new_c_img.stride
    print h, w, s
    shape = (h, s)
    c1 = numpy.ctypeslib.as_array(new_c_img.c1, shape=shape)
    c1[:, :w] = img[:,:,0]

    c2 = numpy.ctypeslib.as_array(new_c_img.c2, shape=shape)
    c2[:, :w] = img[:,:,1]

    c3 = numpy.ctypeslib.as_array(new_c_img.c3, shape=shape)
    c3[:, :w] = img[:,:,2]


def calc_flow(img0, img1):
    h, w, c = img0.shape
    wx = lib.image_new(w, h)
    wy = lib.image_new(w, h)

    im1 = lib.color_image_new(w, h)
    im2 = lib.color_image_new(w, h)

    fill_colot_image_t(img0, im1.contents)
    fill_colot_image_t(img1, im2.contents)

    match_x, match_y, match_z = [ctypes.POINTER(image_t)(),] *3

    params = optical_flow_params_t()
    lib.optical_flow_params_default(ctypes.byref(params))

    lib.optical_flow(wx, wy, im1, im2, params, match_x, match_y, match_z);


    u_computed = numpy_from_image_t(wx.contents)
    v_computed = numpy_from_image_t(wy.contents)

    lib.image_delete(wx)
    lib.image_delete(wy);
    lib.image_delete(match_x);
    lib.image_delete(match_y);
    lib.image_delete(match_z);
    lib.color_image_delete(im1);
    lib.color_image_delete(im2);

    return u_computed, v_computed
