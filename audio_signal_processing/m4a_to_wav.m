source_path = "C:\Users\lenovo\Desktop\Hofer����ת��\1800.m4a";
dest_path = "C:\Users\lenovo\Desktop\Hofer����ת��\wavtest.wav";
[data,fs] = audioread(source_path);
audiowrite(dest_path,data,fs);
