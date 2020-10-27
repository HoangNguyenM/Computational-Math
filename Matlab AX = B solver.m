clc;
clear all;

A = importdata('matrix.csv');
B = importdata('d.csv');

X = linsolve(A,B);

csvwrite('s.csv', X);
