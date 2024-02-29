function [] = display_FNC(sFNC, lw, mx)
%function [] = display_FNC(sFNC, mx)
%sFNC = 53x53 FNC matrix
%mx = scale the data as:  [-mx  +mx], if not present ..
% will set mx to .8 of data max

%%
FT  = readtable('MatchTable_High_NetworkLabling_20190301_ZN.xlsx', 'Sheet','Sheet1', 'Range','A1:K101');
ICN_idx = 10;
ICN_SC = FT{strcmp(FT{:,ICN_idx},'SCN'),2};

ICN_AD = FT{strcmp(FT{:,ICN_idx},'AUD'),2};

ICN_SM = FT{strcmp(FT{:,ICN_idx},'SMN'),2};

ICN_VS = FT{strcmp(FT{:,ICN_idx},'VIS'),2};

ICN_CC = FT{strcmp(FT{:,ICN_idx},'CON'),2};

ICN_DM = FT{strcmp(FT{:,ICN_idx},'DMN'),2};

ICN_CB = FT{strcmp(FT{:,ICN_idx},'CER'),2};

temp_NAN = find(strcmp(FT{:,ICN_idx},'NAN'));

select_ICN = [ICN_SC; ICN_AD; ICN_SM; ICN_VS; ICN_CC; ICN_DM; ICN_CB];
domain_ICN  = {ICN_SC, ICN_AD, ICN_SM, ICN_VS, ICN_CC, ICN_DM, ICN_CB};
num_ICN  = length(select_ICN);
domain_Name = {'SC', 'AUD', 'SM', 'VS', 'CC', 'DM', 'CB'};

%% draw

if (~exist('mx','var'))
    mx = max(abs(sFNC-diag(diag(sFNC))),[],'all','omitmissing');
end

Draw_FNC_Trends(sFNC, domain_Name, domain_ICN, lw,  [-mx mx])

axis image
