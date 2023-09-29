A = readmatrix("Infected data_1.csv"); % reading in the data
Prevalencedata= A(:,1); % first column

Fitted_Parameters=[0.0001 0.2 0.03]; %beta, gamma, alpha %Peak at day 109.
Initial_Guess=Fitted_Parameters; %beta, gamma, alpha. 

%Bounds for the parameters
Lowerbounds = [0 0 0];
Upperbounds=[1 1 1];

options=optimset('Disp','off','TolX',1e-8,'TolFun',1e-8,'MaxIter',15000,'MaxFunEval',15000);  
[EstimatedParameters,fval,exitflag]=fmincon(@(Initial_Guess)err_in_dataSEIR_GLS(Initial_Guess,Prevalencedata),Initial_Guess,[],[],[],[],Lowerbounds, Upperbounds, [], options);