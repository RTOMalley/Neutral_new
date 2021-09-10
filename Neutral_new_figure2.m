%
%   program to generate population using neutral theory
%

% -----------------------------------------------------------------------
%
%   Figure 2 shows results for one population size
%      M = 10,000
%
%   which has 200 different species
%      N = 200
%
%   there is an implied number of individuals per species at generation 0
%      indiviuduals per species = 50
%
%   the 200 species are evenly distributed over ten growth rates
%
%   the growth rates are exponentially distributed from ~0.6 to 1.0
%
%   there is also a rate at which new individuals are brought in
%      s_rate = number brought in every generation
%
%   in the paper this is shown as a percent of the total population
%      percent ~ 0.3
%
%   be sure to change the file name run # to repeat model output
%
% -----------------------------------------------------------------------

clear all

%   numgen = number of generations you want this to run
%   M = total population (constant)
%   N = number of different species
%   P(j) = growth rate array
%   s_rate = number of individuals brought in every generation
%   S = species array
%   SS = species array for reproduction

% ----------------------------
%   open csv file for output
% ----------------------------

fid = fopen ('species200_years10_srate32_run1.csv','w');

numgen=365*10;  %number of generations to run model
M=10000;        %total population
N=200;          %number of different species

% ------------------------------------------------
%   create exponential separation in growth rate
%   for the 200 species
% ------------------------------------------------

gmax = 1.0;     %max growth rate

% ------------------------------------------
%   make 20 species for each division rate
% ------------------------------------------

P(N) = gmax;
for k = 1:19
    P(N-k) = P(N);
end

delta = 0.001;  %  initial delta to separate growth rates

for n=N-20:-20:20
   P(n)= P(n+1)-delta;

   for k = 1:19
     P(n-k) = P(n);
   end

   delta = delta / 0.518;   % modify for next set of growth rates

end

% -------------------------------------
%   introduce speciation rate
%   32 of 10,000 ~ 0.3% of population
% -------------------------------------

s_rate = 32;   %  s_rate = number of new individuals brought in every generation

% -----------------------------
%   initialize species arrays
% -----------------------------

S=ones(1,M);
for i=1:N
    S([(i-1)*M/N+1:i*M/N])=i;
end

SS=zeros(1,2*M);

% -----------------------------------
%  output header and starting point
%  for generation zero
% -----------------------------------

spec0 = unique(S);

fprintf(fid,'generation,species count,');
fprintf(fid,'S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,');
fprintf(fid,'S11,S12,S13,S14,S15,S16,S17,S18,S19,S20,');
fprintf(fid,'S21,S22,S23,S24,S25,S26,S27,S28,S29,S30,');
fprintf(fid,'S31,S32,S33,S34,S35,S36,S37,S38,S39,S40,');
fprintf(fid,'S41,S42,S43,S44,S45,S46,S47,S48,S49,S50,');
fprintf(fid,'S51,S52,S53,S54,S55,S56,S57,S58,S59,S60,');
fprintf(fid,'S61,S62,S63,S64,S65,S66,S67,S68,S69,S70,');
fprintf(fid,'S71,S72,S73,S74,S75,S76,S77,S78,S79,S80,');
fprintf(fid,'S81,S82,S83,S84,S85,S86,S87,S88,S89,S90,');
fprintf(fid,'S91,S92,S93,S94,S95,S96,S97,S98,S99,S100,');
fprintf(fid,'S101,S102,S103,S104,S105,S106,S107,S108,S109,S110,');
fprintf(fid,'S111,S112,S113,S114,S115,S116,S117,S118,S119,S120,');
fprintf(fid,'S121,S122,S123,S124,S125,S126,S127,S128,S129,S130,');
fprintf(fid,'S131,S132,S133,S134,S135,S136,S137,S138,S139,S140,');
fprintf(fid,'S141,S142,S143,S144,S145,S146,S147,S148,S149,S150,');
fprintf(fid,'S151,S152,S153,S154,S155,S156,S157,S158,S159,S160,');
fprintf(fid,'S161,S162,S163,S164,S165,S166,S167,S168,S169,S170,');
fprintf(fid,'S171,S172,S173,S174,S175,S176,S177,S178,S179,S180,');
fprintf(fid,'S181,S182,S183,S184,S185,S186,S187,S188,S189,S190,');
fprintf(fid,'S191,S192,S193,S194,S195,S196,S197,S198,S199,S200,\n');

fprintf(fid,'%3d,%5d',0,length(spec0));

for s = 1:N
    valper = find(S==s);
    count = length(valper);
    fprintf(fid,',%5d',count);
end
fprintf(fid,'\n');

% -----------------------------------
%   iterate over numgen generations
% -----------------------------------

for i=1:numgen

    i  % indicate generation to user

    % ----------------------------------------
    %   let the current population reproduce
    % ----------------------------------------

    SS = S;
    left = M;
    for j=1:N                   % for each species
        I=find(S==j);           % reproduce based on species count and growth rate.
        numtot=P(j)*length(I);  % get the number of their offspring
        num = round(numtot);
        for n = 1:num           % for each offspring
           SS(left+n)=j;
        end
        left = left + num;
    end

    % ------------------------------
    %   choose M survivors from SS
    % ------------------------------

    x = randsample(left,M);  % choose the M that will live from SS
    S(1:M)=SS(x);            % set up surviving generation

    % ---------------------------------------
    %   select new species to import
    % 
    %   use current population distribution
    %   to select species to import 
    %   from neighbors
    % ---------------------------------------

    for g=1:s_rate
      val1 = randsample(M,1);  %  choose one new growth rate
      species = S(val1);
      rate = P(species);
      iset = find(rate==P);    % find original species subset w/ that rate
      val2 = randsample(length(iset),1);  % randomly choose one of them for import
      import = iset(val2);
      val3 = randsample(M,1);  %  choose one to go away
      S(val3) = import;        %  replace it with new species number
    end

    SS=zeros(1,2*M);

    % ----------------------------------------------------
    %   track status of all species over all generations
    % ----------------------------------------------------

    for j=1:N
        I=find(S==j);
        B(i,j)=length(I);
        clear I;
    end
    mu(i)=sum(B(i,:).*P)/sum(B(i,:));
    tester = unique(S);
    numspec(i) = length(tester);

    % ---------------------------------
    %   output the number of species
    %   for this generation
    % ---------------------------------

    spec0 = unique(S);

    fprintf(fid,'%5d,%4d',i,length(spec0));

    for s = 1:N
        valper = find(S==s);
        count = length(valper);
        fprintf(fid,',%5d',count);
    end
    fprintf(fid,'\n');

end

% -----------------------------
%   output simple figure 
%   to see how things changed
%   over the generations
% -----------------------------

figure 
subplot(3,1,1)
plot(1:max(size(B)),B)
ylabel count
subplot(3,1,2)
plot(1:max(size(B)),mu)
ylabel 'avg growth rate'
subplot(3,1,3)
plot(1:max(size(B)),numspec)
ylabel species

fclose(fid)

