WRITE "SIR_PATCH"$

B_:={y1,y2,S1,I1,R1,S2,I2,R2}$
FOR EACH EL_ IN B_ DO DEPEND EL_,T$

B1_:={p11,p22,N1,N2,beta1,beta2,gamma}$

%Number states
NX_:=6$
%Number INPUTS
NU_:=0$
%NUMBER OUTPUTS
NY_:=2$

%FIXING STUFF
LET N1=1000$
LET N2=1000$

%MODEL EQUATIONS
C_:={df(S1,t)=-beta1*p11*S1*(p11*I1/N1+p21*I2/N2)-beta2*p12*S1*(p12*I1/N1+p22*I2/N2),
df(I1,t)=beta1*p11*S1*(p11*I1/N1+p21*I2/N2)+beta2*p12*S1*(p12*I1/N1+p22*I2/N2)-gamma*I1,
df(S2,t)=-beta1*p21*S2*(p11*I1/N1+p21*I2/N2)-beta2*p22*S2*(p12*I1/N1+p22*I2/N2),
df(I2,t)=beta1*p21*S2*(p11*I1/N1+p21*I2/N2)+beta2*p22*S2*(p12*I1/N1+p22*I2/N2)-gamma*I2,
df(R1,t)=gamma*I1,
df(R2,t)=gamma*I2,
y1=R1,
y2=R2
}$

SEED_:=55$
DAISY()$

END$