#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 23 12:54:22 2022

@author: manugv
"""


class DataCont:
    pass


class Dict2Class:
    def __init__(self, arg_dict):
        self.__dict__.update(arg_dict)
