function [ x ] = p_i( n, i )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
x = factorial(n)/(factorial(i)*factorial(n-i));
end

