%ʵ��Ҫ�������ͬ����������ʾ
clc
clear all
close all
N=32;nn=0:(N-1);
subplot(311);
w = ones(N,1);                                          %���δ�ʵ��
stem(nn,w)
xlabel('����');ylabel('����');title('(a)���δ�')
subplot(312);
 w = 0.54 - 0.46*cos(2*pi*(0:N-1)'/(N-1));     %������ʵ��
stem(nn,w)
xlabel('����');ylabel('����');title('(b)������')
subplot(313)
w = 0.5*(1 - cos(2*pi*(0:N-1)'/(N-1)));     %������ʵ��
stem(nn,w)
xlabel('����');ylabel('����');title('(c)������')
