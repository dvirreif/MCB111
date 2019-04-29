%Mi-seq data loaded, those that are commented out were replicates not used in our data
%analysis

AID_minus = fastqread('AID-minus-1-AR.fastq');
%AR_AID_minus2 = fastqread('AID-minus-2-AR.fastq');
AID_plus = fastqread('AID-plus-1-AR.fastq');
%AR_AID_plus2 = fastqread('AID-plus-2-AR.fastq');
BEd_minus = fastqread('BEd-minus-1-AR.fastq');
%AR_BEd_minus2 = fastqread('BEd-minus-2-AR.fastq');
BEd_plus = fastqread('BEd-plus-1-AR.fastq');
%AR_BEd_plus2 = fastqread('BEd-plus-2-AR.fastq');
BEN_minus = fastqread('BEN-minus-1-AR.fastq');
%AR_BEN_minus2 = fastqread('BEN-minus-2-AR.fastq');
BEN_plus = fastqread('BEN-plus-1-AR.fastq');
%AR_BEN_plus2 = fastqread('BEN-plus-2-AR.fastq');
Cas9_minus = fastqread('Cas9-minus-1-AR.fastq');
%AR_Cas9_minus2 = fastqread('Cas9-minus-2-AR.fastq');
Cas9_plus = fastqread('Cas9-plus-1-AR.fastq');
%AR_Cas9_plus2 = fastqread('Cas9-plus-2-AR.fastq');
CDA_minus = fastqread('CDA-minus-1-AR.fastq');
%AR_CDA_minus2 = fastqread('CDA-minus-2-AR.fastq');
CDA_plus = fastqread('CDA-plus-1-AR.fastq');
%AR_CDA_plus2 = fastqread('CDA-plus-2-AR.fastq');
eA3A_minus = fastqread('eA3A-minus-1-AR.fastq');
%AR_eA3A_minus2 = fastqread('eA3A-minus-2-AR.fastq');
eA3A_plus = fastqread('eA3A-plus-1-AR.fastq');
%AR_eA3A_plus2 = fastqread('eA3A-plus-2-AR.fastq');
evo3_minus = fastqread('evo3-minus-1-AR.fastq');
%AR_evo3_minus2 = fastqread('evo3-minus-2-AR.fastq');
%AR_evo3_plus1 = fastqread('evo3-plus-1-AR.fastq');
evo3_plus = fastqread('evo3-plus-2-AR.fastq');
%NO_minus1 = fastqread('NO-minus-1-AR.fastq');
NO_minus = fastqread('NO-minus-2-AR.fastq');
NO_plus = fastqread('NO-plus-1-AR.fastq');
%AR_NO_plus2 = fastqread('NO-plus-2-AR.fastq')

%% Get all data organized into usable format
all_vars = string(who);
len_allvars = length(all_vars);
A = "holder";
for i =1 :len_allvars
    temp = eval(all_vars(i));
    temp = struct2cell(temp);
    temp = string(temp(2,:,:));
    v = genvarname('A',who);
    eval([v ' = squeeze(temp)'])
end

%% Repeats looking for
Motif_5repeats = ['CAGCAGCAGCAGCAG'];
Motif_8repeats = ['CAGCAGCAGCAGCAGCAGCAGCAG']; % arbitratily healthy length as relative repeat length

%% Finding number of repeats in each read

AR_distr = zeros(50,17);
for j = 1:16
    currvar = eval("A" + j);
    len_currvar = length(currvar);
    numrepeats = 8;
    repeats_ones = ones(len_currvar,1);
    motif_repeats = Motif_8repeats;
    for i = 1 :50 
        var_contains = contains(currvar, motif_repeats);
        repeat_ones = repeats_ones.*var_contains;
        num_currentrepeat = sum(repeat_ones);
        AR_distr(i,1) = numrepeats;
        AR_distr(i,j+1) = num_currentrepeat;        
        if num_currentrepeat == 0
            break
        end
        motif_repeats = motif_repeats + "CAG";
        numrepeats = numrepeats + 1;

    end
end

%%%%%%%%%!!!!!!!!!!!!!!!!!!!!! We need to subtract each number of reads with
%%%%%%%%%whatever has more reads, because it is included in all of the
%%%%%%%%%counts (for instance: 5-repeats sequence will be counted from every 
%%%%%%%%%8-repeats sequence etc


AR_real_counts=zeros(50,17);
for i=2:50
    for j = 1:17
        AR_real = (AR_distr(i-1,j))-(AR_distr(i,j));
        AR_real_counts(i-1,j)= AR_real;
        AR_real_counts(i-1,1)= AR_distr(i-1,1) ;
    end
end



%we need to normalize each set of data to the total number of reads for
%each fasta file, becuase number of reads vary. 
%%matrix with all sums. but to be honest, we need to normalize to reads
%%which contain the sequence of the gene of interest, which will be our
%%total

Motif_ARgene = ['GCGAAGTGATCCAGAACCCGGGCCCCAGGCACCCAGAGGCCGCGAGCGCAGCACCTCCCGGCGCCAGTTTGCTGCTGCTG'];

AR_geneseq = zeros(1,17);
for j = 1:16
    currvar = eval("A" + j);
    len_currvar1 = length(currvar);
    repeats_ones1 = ones(len_currvar1,1);
    motif_repeats1 = Motif_ARgene;
    for i = 1 :50 
        var_contains1 = contains(currvar, motif_repeats1);
        repeat_ones1 = repeats_ones1.*var_contains1;
        num_currentrepeat1 = sum(repeat_ones1);
        AR_geneseq(1,j+1) = num_currentrepeat1;
        AR_geneseq(1,1) = 0;
    end
end


AR_distr_norm = zeros(50,17);
for i = 1 :50 
    for j = 1:17
        AR_ratio= (AR_real_counts(i,j))./(AR_geneseq(1,j));
        AR_distr_norm(i,j) = AR_ratio;
        AR_distr_norm(i,1) = AR_distr(i,1);
    end
end



%We need to do some normalization to sgRNA-. I subtracted sgRNA- fractions 
%of given repeats (in theory: all repeats existing) with sgRNA+ results for each BE. It gves us
%delta_counts_of_repeats. The more positive is the difference - the bigger
%is the editing (fraction edited).

AR_distr_norm2 = zeros(50,17);
for i = 1 :50 
    for j = [2 4 6 8 10 12 14 16]
    AR_sgRNA_sub = (AR_distr_norm(i,j))-(AR_distr_norm(i,j+1));
    AR_distr_norm2(i,j) = AR_sgRNA_sub;
    AR_distr_norm2(i,1)=AR_distr(i,1);
    end
end

% new array with repeats and substrated ratio only
% numrepeats AID BEN BEd CDA Cas9 NO eA3A evo3 

AR_distr_norm_changes=(AR_distr_norm2(:, [1 2 4 6 8 10 12 14 16]));
%% Make figures comparing each normalized base editor to no BE control

%we should plot distributions for each given base editor together with NO
%base editor data
figure
scatter(AR_distr_norm(:,1), AR_distr_norm_changes(:,2), 'g', 'filled'); %AID
hold on;
scatter(AR_distr_norm(:,1), AR_distr_norm_changes(:,7), 'r'); %NO
hold off
xlim([5 50])
title('Editing of CAG repeats by chosen base editors')
xlabel('number of CAG repeats')
ylabel('fraction of reads where editing of the target locus occured')
legend({'y = AID','y = no BE'})

figure
scatter(AR_distr_norm(:,1), AR_distr_norm_changes(:,3), 'g', 'filled'); %BEN
hold on;
scatter(AR_distr_norm(:,1), AR_distr_norm_changes(:,7), 'r'); %NO
hold off
xlim([5 50])
title('Editing of CAG repeats by chosen base editors')
xlabel('number of CAG repeats')
ylabel('fraction of reads where editing of the target locus occured')
legend({'y = BEnuclease','y = no BE'})

figure
scatter(AR_distr_norm(:,1), AR_distr_norm_changes(:,4), 'g', 'filled'); %BEd
hold on;
scatter(AR_distr_norm(:,1), AR_distr_norm_changes(:,7), 'r'); %NO
hold off
xlim([5 50])
title('Editing of CAG repeats by chosen base editors')
xlabel('number of CAG repeats')
ylabel('fraction of reads where editing of the target locus occured')
legend({'y = BEdead','y = no BE'})
%% Make figures comparing each normalized base editor to no BE control

figure

scatter(AR_distr_norm(:,1), AR_distr_norm_changes(:,5), 'g', 'filled'); %CDA
hold on;
scatter(AR_distr_norm(:,1), AR_distr_norm_changes(:,7), 'r'); %NO
hold off
xlim([5 50])
grid on
title('Editing of CAG repeats by chosen base editors')
xlabel('number of CAG repeats')
ylabel('fraction of reads where editing of the target locus occured')
legend({'y = CDA','y = no BE'})

figure

scatter(AR_distr_norm(:,1), AR_distr_norm_changes(:,6), 'g', 'filled'); %Cas9
hold on;
scatter(AR_distr_norm(:,1), AR_distr_norm_changes(:,7), 'r'); %NO
hold off
xlim([5 50])
grid on
title('Editing of CAG repeats by chosen base editors')
xlabel('number of CAG repeats')
ylabel('fraction of reads where editing of the target locus occured')
legend({'y = Cas9','y = no BE'})

%% Make figures comparing each normalized base editor to no BE control

figure
scatter(AR_distr_norm(:,1), AR_distr_norm_changes(:,8), 'g', 'filled'); %eA3A
hold on;
scatter(AR_distr_norm(:,1), AR_distr_norm_changes(:,7), 'r'); %NO
hold off
xlim([5 50])
title('Editing of CAG repeats by chosen base editors')
xlabel('number of CAG repeats')
ylabel('fraction of reads where editing of the target locus occured')
legend({'y = eA3A','y = no BE'})

figure
scatter(AR_distr_norm(:,1), AR_distr_norm_changes(:,9), 'g', 'filled'); %evo3
hold on;
scatter(AR_distr_norm(:,1), AR_distr_norm_changes(:,7), 'r'); %NO
hold off
xlim([5 50])
title('Editing of CAG repeats by chosen base editors')
xlabel('number of CAG repeats')
ylabel('fraction of reads where editing of the target locus occured')
legend({'y = evo3','y = no BE'})



%%% We could also plot total number of repeats (let's say 20+) for each
%%% base editor/no base editor. To give an idea how many repeats are there.
%%% Here we actually could normalize to no sgRNA conditions.


AR_distr_20plus=zeros(1,16);
for j = 2:17
    AR_sum20 = sum(AR_real_counts(13:end,j));
    AR_distr_20plus(1,j)= AR_sum20;
end

AR_distr_20plus_fraction = AR_distr_20plus./AR_geneseq;


AR_distr_20plus_fraction_norm = zeros(1,17);
for j = [2 4 6 8 10 12 14 16]
    AR_20plus_sgRNA_sub = (AR_distr_20plus_fraction(1,j))-(AR_distr_20plus_fraction(1,j+1));
    AR_distr_20plus_fraction_norm(1,j) = AR_20plus_sgRNA_sub;
    AR_distr_20plus_fraction_norm(1,1)=0;
end

AR_distr_20plus_fraction_changes=(AR_distr_20plus_fraction_norm(1, [2 4 6 8 10 12 14 16]));

c = categorical({'AID','BEN','BEd', 'CDA', 'Cas9', 'no BE', 'eA3A', 'evo3'});
bar(c,AR_distr_20plus_fraction_changes)
grid on
title('Editing of CAG repeats by chosen base editors')
xlabel('Base editor')
ylabel('fraction of reads edited (for 20+ CAG repeat regions)')





