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
lib.writeFlowFile.restype = ctypes.c_int

lib.color_image_load.restype = ctypes.POINTER(color_image_t)
lib.color_image_load.argtypes = [ctypes.c_char_p]
lib.image_delete.restype = None

lib.color_image_new.restype = ctypes.POINTER(color_image_t)
lib.image_new.restype = ctypes.POINTER(image_t)


def create_params(preset="default"):
    """
    Create instance of optical flow parameters according to one of predefined presets
    :param preset: string, available values: default, sintel, middlebury, kitti
    :return:
    """

    dict_generators = dict(default=lib.optical_flow_params_default,
                           sintel=lib.optical_flow_params_sintel,
                           middlebury=lib.optical_flow_params_middlebury,
                           kitti=lib.optical_flow_params_kitti)
    params = optical_flow_params_t()
    func = dict_generators[preset]

    func(ctypes.byref(params))
    return params


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


def numpy_from_lib_image(img_t):
     if type(img_t) == image_t:
          return numpy_from_image_t(img_t)
     if type(img_t) == color_image_t:
          return numpy_from_color_image_t(img_t)
     raise Exception("Unsupported input type")


def read_flow(path):
    t = lib.readFlowFile(path)

    u = numpy_from_image_t(t[0].contents)
    v = numpy_from_image_t(t[1].contents)

    lib.image_delete(t[0])
    lib.image_delete(t[1])

    return u, v


def write_flow(path, u, v):
    """
    Save optical flow (u,v) to .flo file.
    :param path: destination path
    :param u: x component of optical flow
    :param v: y component of optical flow
    :return: None
    """
    if u.shape != v.shape:
        raise Exception("Shapes of x- and y-components of optical flow must match")

    if u.ndim != 2:
        raise Exception("Input optical flow components must me 2-dimensional")

    h, w = u.shape

    wx = lib.image_new(w, h)
    wy = lib.image_new(w, h)

    fill_image_t(u, wx.contents)
    fill_image_t(v, wy.contents)

    res = lib.writeFlowFile(path, wx, wy)

    lib.image_delete(wx)
    lib.image_delete(wy)

    if res != 0:
        raise Exception("Failed to save optical flow")


def fill_colot_image_t(img, new_c_img):
    h, w, s = new_c_img.height, new_c_img.width, new_c_img.stride
    shape = (h, s)
    c1 = numpy.ctypeslib.as_array(new_c_img.c1, shape=shape)
    c1[:, :w] = img[:,:,0]

    c2 = numpy.ctypeslib.as_array(new_c_img.c2, shape=shape)
    c2[:, :w] = img[:,:,1]

    c3 = numpy.ctypeslib.as_array(new_c_img.c3, shape=shape)
    c3[:, :w] = img[:,:,2]


def fill_image_t(img_src, img_dst):
    h, w, s = img_dst.height, img_dst.width, img_dst.stride
    shape = (h, s)
    data = numpy.ctypeslib.as_array(img_dst.data, shape=shape)
    data[:, :w] = img_src[:,:]


def calc_flow(img0, img1, params=None):
    """
    Calculate optical flow between two images. For each pixel in the source image the algorithm searches for
    corresponding location in the reference image.
    :param img0: Source image
    :param img1: Reference image
    :param params: parameters of optical flow algorithm
    :return: (u, v) tuple for horizontal (u) and vertical (v) components of optical flow
    """
    h, w, c = img0.shape
    wx = lib.image_new(w, h)
    wy = lib.image_new(w, h)

    im1 = lib.color_image_new(w, h)
    im2 = lib.color_image_new(w, h)

    fill_colot_image_t(img0, im1.contents)
    fill_colot_image_t(img1, im2.contents)

    match_x, match_y, match_z = [ctypes.POINTER(image_t)(),] * 3

    if params is None:
        params = optical_flow_params_t()
        lib.optical_flow_params_default(ctypes.byref(params))

    lib.optical_flow(wx, wy, im1, im2, ctypes.byref(params), match_x, match_y, match_z)

    u_computed = numpy_from_image_t(wx.contents)
    v_computed = numpy_from_image_t(wy.contents)

    lib.image_delete(wx)
    lib.image_delete(wy)
    lib.image_delete(match_x)
    lib.image_delete(match_y)
    lib.image_delete(match_z)
    lib.color_image_delete(im1)
    lib.color_image_delete(im2)

    return u_computed, v_computed


def warp_image(image, u, v):
    shape = list(image.shape)
    if len(shape) == 2:
        shape.append(1)
    elif len(shape) != 3:
         raise Exception("Unsupported number of dimensions")
    h, w, channels = shape

    if channels == 3:
         func_new_image = lib.color_image_new
         func_fill_image = fill_colot_image_t
         func_delete_image = lib.color_image_delete
         func_warp = lib.color_image_warp
    elif channels == 1:
         func_new_image = lib.image_new
         func_fill_image = fill_image_t
         func_delete_image = lib.image_delete
         func_warp = lib.image_warp
    else:
         raise Exception("Unsupported format")
         
    
    im_warped = func_new_image(w, h)
    im        = func_new_image(w, h)
    
    mask = lib.image_new(w, h)
    wx = lib.image_new(w, h)
    wy = lib.image_new(w, h)
    
    fill_image_t(u, wx.contents)
    fill_image_t(v, wy.contents)
    
    func_fill_image(image, im.contents)
    
    func_warp(im_warped, mask, im, wx, wy)
    im_warped_np = numpy_from_lib_image(im_warped.contents)
    
    lib.image_delete(wx)
    lib.image_delete(wy);
    lib.image_delete(mask);

    func_delete_image(im_warped);
    func_delete_image(im);

    return im_warped_np
    
