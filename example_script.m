% script file to illustrate the determination of a trend in the MgII index



load MgII_index_data
clf
clc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n')
disp('******** First plot raw MgII data ********')
disp('-- Hit any key to continue --')
pause

clf
h = plot(time_MgII,MgII_Bremen,time_MgII,MgII_LASP);
set(h(:),'linewidth',0.8); 
datetick('x','yyyy')
grid on
xlabel('year')
ylabel('MgII index')
legend('Bremen','LASP','Location','northeast')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n')
disp('******** Run the model in manual mode, with order = 12 ********')
disp('-- In the following : the number of levels is 100 -- ')
disp('-- In the following : variations on scales < 81 days are discarded -- ')
disp('-- Hit any key to continue --')
pause

width = 81;
number_levels = 100;
order = 12;
verbose = 1;


[f,df,order,condition,BIC] = solar_trend_estimation(MgII_Bremen,MgII_LASP,width,order,number_levels,verbose);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n')
disp('******** Run the model again, with verbose set to maximum ********')
disp('-- Hit any key to continue --')
pause

width = 81;
number_levels = 100;
order = 12;
verbose = 2;


[f,df,order,condition,BIC] = solar_trend_estimation(MgII_Bremen,MgII_LASP,width,order,number_levels,verbose);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n')
disp('******** Run the model again, with verbose set to maximum ********')
disp('-- The number of levels is set here to 500 : this slows down the --')
disp('-- processing. However, the effective nr of pairs does not change much --')
disp('-- Hit any key to continue --')
pause

width = 81;
number_levels = 500;
order = 12;
verbose = 2;


[f,df,order,condition,BIC] = solar_trend_estimation(MgII_Bremen,MgII_LASP,width,order,number_levels,verbose);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n')
disp('******** Run the model in semi-automatic mode and let the user decide which order is the best ********')
disp('-- Hit any key to continue --')
pause

width = 81;
number_levels = 100;
order = -1;
verbose = 1;

[f,df,order,condition,BIC] = solar_trend_estimation(MgII_Bremen,MgII_LASP,width,order,number_levels,verbose);






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n')
disp('******** Run again the model in fully automatic mode ********')
disp('-- This mode is NOT recommended : it is always better to inspect the data ')
disp('-- Hit any key to continue --')
pause


width = 81;
number_levels = 100;
order = -2;
verbose = 1;

[f,df,order,condition,BIC] = solar_trend_estimation(MgII_Bremen,MgII_LASP,width,order,number_levels,verbose);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n')
disp('******** Show the corrected values of the MgII index from LASP ********')
disp('-- Hit any key to continue --')
pause


[f,df,order,condition,BIC] = solar_trend_estimation(MgII_Bremen,MgII_LASP,width,12,number_levels,0);
MgII_LASP_corrected = MgII_LASP./f;

clf
subplot(211)
h = plot(time_MgII,MgII_Bremen,time_MgII,MgII_LASP);
set(h(:),'linewidth',0.8); 
datetick('x','yyyy')
grid on
xlabel('year')
ylabel('MgII index')
legend('Bremen','LASP-uncorrected','Location','northeast')

subplot(212)
h = plot(time_MgII,MgII_Bremen,time_MgII,MgII_LASP_corrected);
set(h(:),'linewidth',0.8); 
datetick('x','yyyy')
grid on
xlabel('year')
ylabel('MgII index')
legend('Bremen','LASP-corrected','Location','northeast')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n')
disp('******** Show the ratio between the two MgII indices ********')
disp('-- Hit any key to continue --')
pause


[f,df,order,condition,BIC] = solar_trend_estimation(MgII_Bremen,MgII_LASP,width,12,number_levels,0);
MgII_LASP_corrected = MgII_LASP./f;

clf
subplot(211)
h = plot(time_MgII,MgII_LASP./MgII_Bremen,time_MgII,f);
ax = axis;
set(h(:),'linewidth',0.8); 
datetick('x','yyyy')
grid on
xlabel('year')
ylabel('ratio uncorrected LASP / Bremen')
legend('ratio uncorrected LASP / Bremen','correction f','Location','Northwest')

subplot(212)
h = plot(time_MgII,MgII_LASP_corrected./MgII_Bremen);
axis(ax)
set(h(:),'linewidth',0.8); 
datetick('x','yyyy')
grid on
xlabel('year')
ylabel('ratio corrected LASP / Bremen')








