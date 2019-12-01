clear all; format long; clc; close all; %#ok<CLALL>
warning off

d = daq.getDevices;
s = daq.createSession('ni');
s.Rate = 51200;%2^10;%1024;%51200;%51200;
Fs=s.Rate;
NFFT=1*Fs;%2^nextpow2(4*Fs);%4*Fs;%2^15;%2^nextpow2(s.Rate);%2^16; 2^nextpow2(fs)
s.DurationInSeconds = 10; %Nfft/Fs;%Nfft/Fs;

[ch,idx] = addAnalogInputChannel(s,'cDAQ1Mod1',[0, 1, 2, 3],'IEPE'); % canal 0, até 3

d.Subsystems.ChannelNames;
%%
sens1= 1.017/1000;% %sensibilidade V/(m/s2) - acc x
calscale1=1/sens1;
%%
sens2= 1.057/1000;% %sensibilidade V/(m/s2)  - acc y
calscale2=1/sens2;
% %%
% sens3=0.9701/1000;% %sensibilidade V/(m/s2) - acc z
% calscale3=1/sens3;
% %%
% sens4=0.1087;% %sensibilidade V/(N) - Força
% calscale4=1/sens4;
%%
trec=s.DurationInSeconds; % tempo de gravação (depende do TR da sala)
nbits=24;

AMP=1;%0;%0.95;% amplitude do sinal de ruído rosa [0-1]

%%

for nn=1:1%numero de medias
%% SOUND CAPTURE

disp('Start of Recording')
[Recording ,t] = s.startForeground;
%%
disp('End');

acc1IN=calscale1.*Recording(:,1);% canal 1 Acc em baixo do Coxim
acc2OUT=calscale2.*Recording(:,2);% canal 2 Acc em cima do Coxim
% accz=calscale3.*Recording(:,3);% canal 3
% forca=calscale4.*Recording(:,4);% canal 4 

audiowrite('acc1IN.wav',Recording(:,1),Fs); % canal 1
audiowrite('acc2OUT.wav',Recording(:,2),Fs); % canal 2
% audiowrite('accz.wav',Recording(:,3),Fs); % canal 3
% audiowrite('forca.wav',Recording(:,4),Fs); % canal 4

figure()
plot(t,Recording(:,1),'r',t,Recording(:,2),'b'); grid on; %hold on;
xlabel('Tempo [s]'); ylabel('Amplitude [V]')
legend('acc1 IN','acc2 OUT')

clear Recording

a=2; b=1;

[Hf_x, freqh]= tfestimate(acc2OUT,acc1IN,hanning(NFFT),NFFT/a,b*NFFT,Fs);('twosided');
[coe_x, freqc]= mscohere(acc2OUT,acc1IN,hanning(NFFT),NFFT/a,b*NFFT,Fs);

figure()
loglog(freqh,(abs(Hf_x)));  grid on; hold on
xlabel('Frequência [Hz]')
ylabel('Função de Transferência H(f) [m/s2/ m/s2]')
title('Hx(f)')
xlim([20 2000])


figure()
semilogx(freqc,coe_x);  grid on; hold on
xlabel('Frequência [Hz]')
ylabel('Coerência [-]')
title('Coe_z(f)')
% xlim([20 500])
% ylim([0.5 1])

end

% weighting='flat'; % 'A'; 'C', 'flat
% bandtype='third'; % oct, third
% filename='acc1IN.wav';
% Transdutor='acc'; % mic, acc - microfone ou acelerômetro
% Ref=1;% m/s^2
% Nfft=Fs;%2^15; % 2^nextpow2(Fs)
% [fb,Pxxb_accx,fe,Pxx_accx]= Level_vslm(sens1, weighting, bandtype, filename, Fs, Ref, NFFT, Transdutor);
% legend('acc1 IN ')
% 
% weighting='flat'; % 'A'; 'C', 'flat
% bandtype='third'; % oct, third
% filename='acc2OUT.wav';
% Transdutor='acc'; % mic, acc - microfone ou acelerômetro
% Ref=1;% m/s^2
% NFFT=Fs;%2^15; % 2^nextpow2(Fs)
% [~,Pxxb_accy,~,Pxx_accy]= Level_vslm(sens2, weighting, bandtype, filename, Fs, Ref, NFFT, Transdutor);
% legend('acc2 OUT ')
% 
% weighting='flat'; % 'A'; 'C', 'flat
% bandtype='third'; % oct, third
% filename='accz.wav';
% Transdutor='acc'; % mic, acc - microfone ou acelerômetro
% Ref=1;% m/s^2
% NFFT=Fs;%2^15; % 2^nextpow2(Fs)
% [~,Pxxb_accz,~,Pxx_accz]=Level_vslm(sens3, weighting, bandtype, filename, Fs, Ref, NFFT, Transdutor);
% legend('acc z ')

pause
%salva figura em pdf - valores 5 e 4 alteram o formato
figure(1)
set(gcf, 'PaperPosition', [0 0 5 4]); %Position plot at left hand corner with width 5 and height 5. 
set(gcf, 'PaperSize', [5 4]); %Set the paper to have width 5 and height 5. 
saveas(gcf, 'graf1', 'pdf') %Save figure
savefig('graf1.fig')
figure(2)
set(gcf, 'PaperPosition', [0 0 5 4]); %Position plot at left hand corner with width 5 and height 5. 
set(gcf, 'PaperSize', [5 4]); %Set the paper to have width 5 and height 5. 
saveas(gcf, 'graf2', 'pdf') %Save figure
savefig('graf2.fig')
figure(3)
set(gcf, 'PaperPosition', [0 0 5 4]); %Position plot at left hand corner with width 5 and height 5. 
set(gcf, 'PaperSize', [5 4]); %Set the paper to have width 5 and height 5. 
saveas(gcf, 'graf3', 'pdf') %Save figure
savefig('graf3.fig')
% figure(4)
% set(gcf, 'PaperPosition', [0 0 5 4]); %Position plot at left hand corner with width 5 and height 5. 
% set(gcf, 'PaperSize', [5 4]); %Set the paper to have width 5 and height 5. 
% saveas(gcf, 'graf4', 'pdf') %Save figure
% savefig('graf4.fig')
% figure(5)
% set(gcf, 'PaperPosition', [0 0 5 4]); %Position plot at left hand corner with width 5 and height 5. 
% set(gcf, 'PaperSize', [5 4]); %Set the paper to have width 5 and height 5. 
% saveas(gcf, 'graf5', 'pdf') %Save figure
% savefig('graf5.fig')
% figure(6)
% set(gcf, 'PaperPosition', [0 0 5 4]); %Position plot at left hand corner with width 5 and height 5. 
% set(gcf, 'PaperSize', [5 4]); %Set the paper to have width 5 and height 5. 
% saveas(gcf, 'graf6', 'pdf') %Save figure
% savefig('graf6.fig')
% figure(7)
% set(gcf, 'PaperPosition', [0 0 5 4]); %Position plot at left hand corner with width 5 and height 5. 
% set(gcf, 'PaperSize', [5 4]); %Set the paper to have width 5 and height 5. 
% saveas(gcf, 'graf7', 'pdf') %Save figure
% savefig('graf7.fig')
% figure(8)
% set(gcf, 'PaperPosition', [0 0 5 4]); %Position plot at left hand corner with width 5 and height 5. 
% set(gcf, 'PaperSize', [5 4]); %Set the paper to have width 5 and height 5. 
% saveas(gcf, 'graf8', 'pdf') %Save figure
% savefig('graf8.fig')
% figure(9)
% set(gcf, 'PaperPosition', [0 0 5 4]); %Position plot at left hand corner with width 5 and height 5. 
% set(gcf, 'PaperSize', [5 4]); %Set the paper to have width 5 and height 5. 
% saveas(gcf, 'graf9', 'pdf') %Save figure
% savefig('graf9.fig')
% figure(10)
% set(gcf, 'PaperPosition', [0 0 5 4]); %Position plot at left hand corner with width 5 and height 5. 
% set(gcf, 'PaperSize', [5 4]); %Set the paper to have width 5 and height 5. 
% saveas(gcf, 'graf10', 'pdf') %Save figure
% savefig('graf10.fig')
