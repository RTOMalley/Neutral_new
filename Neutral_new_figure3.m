%
%   program to generate population using neutral theory
%

% -----------------------------------------------------------------------
%
%   Figure 3 shows results for two different population sizes
%      M = 10,000 and 100,000
%
%   each has the same number of species
%      N = 200
%
%   there is an implied number of individuals per species
%      indiviuduals per species = 50 and 500
%
%   there is also a rate at which new individuals are brought in
%      m_rate = number brought in every generation
%
%      *** NOTE:  m_rate is equivalent to s_rate in code for Figure 1 and 2
%
%   in the paper m_rate is shown as a percent of the total population
%      percent = 0.03
%
%   because of different population sizes 
%   this makes m_rate change for the different percentages
%
%      M = 10,000
%         m_rate = 3 is 0.03% of 10,000
%
%      M = 100,000 
%         m_rate = 30 is 0.03% of 100,000
%
%   there is also a chance of introducing one new species (0.00002% chance)
%      s_rate = 1 is the one new species brought in
%
%   to get such a low percentage chance, we must look (on average)
%   over multiple generations.  This depends on the population size
%
%      M = 10,000
%         genchance = 500
%
%      M = 100,000 
%         genchance = 50
%
%            to clarify, 1 (s_rate) individual in 50 generations
%            is the same as 2 individuals in 100 generations
%
%            selection 2 individuals in 1 generation for M = 100,000
%            is the same as selecting 0.002 in 100
%
%            inceasing the total possible over 100 generations
%            will shift 0.002 in 100 to a 0.00002% chance
%
%
%   make the changes to appropriate variables for the population size 
%   for the run you wish:
%
%      M = 10,000
%         genchance = 500
%         m_rate = 3
%
%      M = 100,000 
%         genchance = 50
%         m_rate = 30
%
%
%   this includes the naming of the csv output file
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

fid = fopen ('pop100000_mrate30_srate1_in50gen
_run1.csv','w');

numgen=20000;  %number of generations to run model
M=100000;      %total population
N=200;         %number of different species

ireport = 100;  %  output the generation info every ireport generations

% ----------------------------------------------------------------
%   use these keys to randomly sample 1 in 50 times (on average)
%   to introduce a 'neutral' or 'better' species
% ----------------------------------------------------------------

genchance = 50;  %  one chance in genchance to introduce new species  

n_key = 14;  %  random key for inclusion of new 'neutral' species
b_key = 38;  %  random key for inclusion of new 'better' species

deltanew = 0.0010;  %  increase to gmax for better growth rate

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

delta = 0.0010;     %  initial delta to separate growth rates

for n=N-20:-20:20
   P(n)= P(n+1)-delta;

   for k = 1:19
     P(n-k) = P(n);
   end

   delta = delta / 0.518;   % modify for next set of growth rates

end

pindex = N;  % index of maximum species ID (new species will increase this)

% ---------------------------------------
%   introduce mixing rate
%   30 of 100,000 = 0.03% of population
% ---------------------------------------

s_rate = 1;   %  s_rate = number of (new) species every time it is applied
m_rate = 30;  %  set mixing rate at 0.03% = 30/100,000

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

fprintf(fid,'generation,species count\n');
fprintf(fid,'%3d,%5d',0,length(spec0));

fprintf(fid,'\n');

% -----------------------------------
%   iterate over numgen generations
% -----------------------------------

for i=1:numgen

    remval = rem(i,ireport);
    if ( remval == 0 )
        i  % indicate generation to user
    end

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

    x = randsample(left,M);  % choose the 10,000 that will live from SS
    S(1:M)=SS(x);            % set up surviving generation

    % ---------------------------------------
    %   do mixing every generation  
    %
    %   use current population distribution
    %   to select species to import 
    %   from neighbors
    % ---------------------------------------

    for g=1:m_rate
      val1 = randsample(M,1);  %  choose one new growth rate
      species = S(val1);
      rate = P(species);
      iset = find(rate==P);    %  find original species subset w/ that rate
      val2 = randsample(length(iset),1);  % randomly choose one of them for import
      import = iset(val2);
      val3 = randsample(M,1);  %  choose one to go away
      S(val3) = import;        %  replace it with new species number
    end

    % -----------------------------------------------------------------------
    %   bring in a 'neutral' species once every 50 generations (on average)
    % -----------------------------------------------------------------------

    n_test = randsample(genchance,1);

    if ( n_test == n_key )

        % ------------------------------------
        %   select new species growth rate 
        %   use current growth distribution
        %   to select
        % -----------------------------------

        for g=1:s_rate

            val1 = randsample(M,1);  %  choose one current growth rate
            species = S(val1);
            rate = P(species);

            pindex = pindex + 1;     %  prepare to keep track of new species
            P(pindex) = rate;
            N = N+1;
            import = N;

            val3 = randsample(M,1);  %  choose one to go away
            S(val3) = import;        %  replace it with new species number
      end

    end % if (n_test == n_key)

    % ----------------------------------------------------------------------
    %   bring in a 'better' species once every 50 generations (on average)
    % ----------------------------------------------------------------------

    b_test = randsample(genchance,1);

    if ( b_test == b_key )

        % -----------------------------------
        %   select new species growth rate 
        %   use increase in max growth rate
        %   for new growth rate
        % -----------------------------------

        gmax = gmax + deltanew;      %  increase gmax for better growth rate

        for g=1:s_rate
            rate = gmax;

            pindex = pindex + 1;     %  prepare to keep track of new species
            P(pindex) = rate;
            N = N+1;
            import = N;

            val3 = randsample(M,1);  %  choose one to go away
            S(val3) = import;        %  replace it with new species number
      end

    end % if (b_test == b_key)

    SS=zeros(1,2*M);

    % ----------------------------------------------------
    %   track status of all species over all generations
    % ----------------------------------------------------

    tester = unique(S);
    numspec(i) = length(tester);

    % ---------------------------------
    %   output the number of species
    %   for this generation
    % ---------------------------------

    spec0 = unique(S);

    fprintf(fid,'%5d,%4d',i,length(spec0));

    fprintf(fid,'\n');

end

figure
plot(numspec)

fclose(fid)

