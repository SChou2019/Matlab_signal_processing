%ʵ��Ҫ����������˲������
clc
clear all
fs=8000;
%�ѽ�ֹƵ��ת�ɻ��ȱ�ʾ
wp=550*2/fs;
ws=500*2/fs;        
rp=3;
rs=80;
Nn=512;
[n,wn]=ellipord(wp,ws,rp,rs);
[b,a]=ellip(n,rp,rs,wn);
freqz(b,a,Nn,fs);
