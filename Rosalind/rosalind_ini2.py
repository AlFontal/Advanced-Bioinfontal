#! /usr/bin/env python
"""
Author: Alejandro Fontal
Student Registration Number: 920110242090
"""
from sys import argv

def sq_hypotenuse(a, b):
    """
    Calculates the square of an hypotenuse given the length of the legs.

    :param a: Integer, length of leg 1.
    :param b: Integer, length of leg 2.
    :return: h_squared, integer, square of the length of the hypotenuse.
    """

    h_squared = b**2 + a**2
    return h_squared

if __name__ == '__main__':

    input_file = open(argv[1], "r")
    nums = input_file.read().strip()
    num_list = nums.split(" ")

    a = int(num_list[0])
    b = int(num_list[1])
    h_squared = sq_hypotenuse(a, b)
    print h_squared