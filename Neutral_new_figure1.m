%
%   program to generate population using neutral theory
%

% -----------------------------------------------------------------------
%
%   Figure 1 shows results for three different population sizes
%      M = 10,000, 100,000 and 1,000,000
%
%   each has the same number of species
%      N = 10,000
%
%   there is an implied number of individuals per species
%      indiviuduals per species = 1, 10, and 100
%
%   there is also a rate at which new individuals are brought in
%      s_rate = number brought in every generation
%
%   in the paper this is shown as a percent of the total population
%      percent = 0, 0.03 and 0.3
%
%   because of different population sizes 
%   this makes s_rate change for the different percentages
%
%      M = 10,000
%         s_rate = 0 is 0% of 10,000
%         s_rate = 3 is 0.03% of 10,000
%         s_rate = 30 is 0.3% of 10,000
%
%      M = 100,000 
%         s_rate = 0 is 0% of 100,000
%         s_rate = 30 is 0.03% of 100,000
%         s_rate = 300 is 0.3% of 100,000
%
%      M = 1,000,000 
%         s_rate = 0 is 0% of 1,000,000
%         s_rate = 300 is 0.03% of 1,000,000
%         s_rate = 3,000 is 0.3% of 1,000,000
%
%   make the changes to appropriate variables below for the run 
%   you are making, including the naming of the csv output file
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

fid = fopen ('species_1e4_indivs_1_srate_3_run1.csv','w');

numgen=100000;  %number of generations to run model
M=10000;        %total population
N=10000;        %number of different species

% -------------------------------------------
%   create same growth rate for all species
% -------------------------------------------

gmax = 1.0;     %max growth rate
P(1:N) = gmax;  %growth rates for N different species

% -------------------------------------
%   introduce speciation rate
%   3 of 10,000 = 0.03% of population
% -------------------------------------

s_rate = 3;   %  s_rate = number of new individuals brought in every generation

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

   fprintf(fid,'generation,species count,s_rate = 3 = 0.03 percent\n');
   fprintf(fid,'%3d,%5d\n',0,length(spec0));

% -----------------------------------
%   iterate over numgen generations
% -----------------------------------

for i=1:numgen 

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

    % ----------------------------------------
    %   select new species to import
    % 
    %   use current population distribution
    %   to select species to import 
    %   from neighbors
    %  ------------------------------------------------
    %   this is set up for a generalized case
    %   but in this case there is only one growth rate
    % --------------------------------------------------

    for g=1:s_rate
        val1 = randsample(M,1);  %  choose one new growth rate
        species = S(val1);
        rate = P(species);
        iset = find(rate==P);    % find original species subset w/ that rate
        val2 = randsample(length(iset),1); % randomly choose one for import
        import = iset(val2);
        val3 = randsample(M,1);  %  choose one to go away
        S(val3) = import;        %  replace it with new species number
    end

    SS=zeros(1,2*M);

    % -------------------------------------------
    %   no need to output every generation
    %
    %   set up to output often at the beginning
    %   then shift to longer periods
    % --------------------------------------------

    if (i <= 30 )
        rtest = 5;
    end
    if ((i > 30)&(i<=100))
        rtest = 10;
    end
    if ((i>100)&(i<=3000))
        rtest = 25;
    end
    if ((i>3000))
        rtest = 100;
    end

    % --------------------------------------
    %   see if it's time to output results
    % --------------------------------------

    remval = rem(i,rtest);
    if ( remval == 0 )

        i  % indicate generation to user

        % ---------------------------------
        %   output the number of species
        %   for this generation
        % ---------------------------------

        spec0 = unique(S);

        fprintf(fid,'%3d,%5d\n',i,length(spec0));

        length(spec0)  % indicate number of remaining species to user

    end

    % -------------------------------------
    %   if there is only one species left
    %   then nothing will change
    % -------------------------------------

    if (length(spec0) == 1)
        break;
    end

% -------------------------------------------------------
%   you may want to activate the matfile save below
%   if you are working with M = 1,000,000
%   (this is a long run)
%
%   the matfile will make it 
%   so you can restart your program
%   if your system hangs and you need to reboot
%
%   if this happens, 
%   you can modify the code above to continue
%   with the current state saved in the matfile 
%   (rather than starting over)
% -------------------------------------------------------
%
%    remval2 = rem(i,5000);
%    if ( remval2 == 0 )
%        save species_1e4_indivs_100_srate_0.mat
%    end

end

fclose(fid)

