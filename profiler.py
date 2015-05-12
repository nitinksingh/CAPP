# -*- coding: utf-8 -*-
"""
A collection of time profileing decorators from Zapier Enggineering blog

Created on Sat May  2 11:42:25 2015

@author: nitin
"""
from __future__ import division
import time
import cProfile
from line_profiler import LineProfiler

# cProfiling decorator
def do_cprofile(func):
    def profiled_func(*args, **kwargs):
        profile = cProfile.Profile()
        try:
            profile.enable()
            result = func(*args, **kwargs)
            profile.disable()
            return result
        finally:
            profile.print_stats()
    return profiled_func


# Simple time decorator
def timefunc(f):
    def f_timer(*args, **kwargs):
        start = time.time()
        result = f(*args, **kwargs)
        end = time.time()
        print f.__name__, 'took', end - start, 'secs'
        return result
    return f_timer


# Do line by line profiling for function specified in the argument
def do_profile(follow=[]):
        def inner(func):
            def profiled_func(*args, **kwargs):
                try:
                    profiler = LineProfiler()
                    profiler.add_function(func)
                    for f in follow:
                        profiler.add_function(f)
                    profiler.enable_by_count()
                    return func(*args, **kwargs)
                finally:
                    profiler.print_stats()
            return profiled_func
        return inner

    

def func2():
    a = [i**3 for i in range(1000)] 

# Example
@timefunc
@do_profile(follow=[func2])
def func1():
    b = [i**2 for i in range(1000)] 
    func2()
    
if __name__=='__main__':
    func1()