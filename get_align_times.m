function [align_times, trial_numbers] = get_align_times(filters, trial_times, trial_words, alignment)

    % get_align_times.m
    %
    % function that returns aligment times for specified event type in MAD sEEG
    % data
    %
    % INPUTS:
    % filters           structure specifying trial parameters (from Setup.mat
    %                   file)
    % trial_times       event times for every trial in TRIAL# X 1 cell array 
    %                   (from Setup.mat file)
    % trial_words       event codes for every trial in TRIAL# X 1 cell array 
    %                   (from Setup.mat file)
    % alignment         string specifying alignment type, e.g. 'response',
    %                   'trialstart', 'win', 'loss', etc.
    %
    % OUTPUTS:
    % align_times       vector of times for alignment when specified event type
    %                   occurred
    % trial_numbers     vector of trial numbers corresponding to times in
    %                   align_times
    %
    % Author: Aaron Sampson <asampson@jhu.edu>
    % Dec 2021;
        
    align_times = [];    
    trial_numbers = [];
    for trial_num = 1:numel(trial_times)
        single_trial_times = trial_times{trial_num};
        trial_codes = trial_words{trial_num};
        trial_number = filters.trial(trial_num);
        switch alignment
            case 'trialstart'
                align_times = [align_times; single_trial_times(trial_codes==1)];
                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==1)))];
            case 'inspection'
                align_times = [align_times; single_trial_times(trial_codes>17&trial_codes<34)];
                trial_numbers = [trial_numbers; trial_number*ones(size(trial_codes>17&trial_codes<34))];
            case 'amount'
                align_times = [align_times; single_trial_times(mod(trial_codes,2)==0&trial_codes>17&trial_codes<34)];
                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(mod(trial_codes,2)==0&trial_codes>17&trial_codes<34)))];
            case 'probability'
                align_times = [align_times; single_trial_times(mod(trial_codes,2)==1&trial_codes>17&trial_codes<34)];
                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(mod(trial_codes,2)==1&trial_codes>17&trial_codes<34)))];
            case 'tap'
                align_times = [align_times; single_trial_times(trial_codes>1&trial_codes<18)];
                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes>1&trial_codes<18)))];
            case 'nonmotor'
                align_times = [align_times; single_trial_times(trial_codes==36)+0.5];
                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==36)))];
            case 'win_domain'
                align_times = [align_times; single_trial_times(trial_codes==18|trial_codes==19|trial_codes==20|trial_codes==21|trial_codes==26|trial_codes==27|trial_codes==30|trial_codes==31)];
                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==18|trial_codes==19|trial_codes==20|trial_codes==21|trial_codes==26|trial_codes==27|trial_codes==30|trial_codes==31)))];
            case 'loss_domain'
                align_times = [align_times; single_trial_times(trial_codes==22|trial_codes==23|trial_codes==24|trial_codes==25|trial_codes==28|trial_codes==29|trial_codes==32|trial_codes==33)];
                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==22|trial_codes==23|trial_codes==24|trial_codes==25|trial_codes==28|trial_codes==29|trial_codes==32|trial_codes==33)))];
            case 'win_domain_trialstart'
                if ~isnan(filters.opt1_win_amount(trial_num))&&isnan(filters.opt1_loss_amount(trial_num))
                    align_times = [align_times; single_trial_times(trial_codes==1)];
                    trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==1)))];
                else
                    continue
                end
            case 'loss_domain_trialstart'
                if isnan(filters.opt1_win_amount(trial_num))&&~isnan(filters.opt1_loss_amount(trial_num))
                    align_times = [align_times; single_trial_times(trial_codes==1)];
                    trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==1)))];
                else
                    continue
                end
            case 'response'
                align_times = [align_times; single_trial_times(trial_codes==36)];
                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==36)))];
            case 'win_domain_response'
                if ~isnan(filters.opt1_win_amount(trial_num))&&isnan(filters.opt1_loss_amount(trial_num))
                    align_times = [align_times; single_trial_times(trial_codes==36)];
                    trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==36)))];
                else
                    continue
                end
            case 'loss_domain_response'
                if isnan(filters.opt1_win_amount(trial_num))&&~isnan(filters.opt1_loss_amount(trial_num))
                    align_times = [align_times; single_trial_times(trial_codes==36)];
                    trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==36)))];
                else
                    continue
                end
            case 'select_top'
                if filters.horiz_vert(trial_num)==1&&filters.choice(trial_num)==1
                    align_times = [align_times; single_trial_times(trial_codes==36)];
                    trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==36)))];
                else
                    continue
                end
            case 'select_bottom'
                if filters.horiz_vert(trial_num)==1&&filters.choice(trial_num)==2
                    align_times = [align_times; single_trial_times(trial_codes==36)];
                    trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==36)))];
                else
                    continue
                end
            case 'select_left'
                if filters.horiz_vert(trial_num)==0&&filters.choice(trial_num)==1
                    align_times = [align_times; single_trial_times(trial_codes==36)];
                    trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==36)))];
                else
                    continue
                end
            case 'select_right'
                if filters.horiz_vert(trial_num)==0&&filters.choice(trial_num)==2
                    align_times = [align_times; single_trial_times(trial_codes==36)];
                    trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==36)))];
                else
                    continue
                end
            case 'outcome'
                align_times = [align_times; single_trial_times(trial_codes==37)];
                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==37)))];
            case 'win'
                align_times = [align_times; single_trial_times(trial_codes==37&filters.result(trial_num)>0)];
                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==37&filters.result(trial_num)>0)))];
            case 'loss'
                align_times = [align_times; single_trial_times(trial_codes==37&filters.result(trial_num)<0)];
                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==37&filters.result(trial_num)<0)))];
            case 'neutral'
                align_times = [align_times; single_trial_times(trial_codes==37&filters.result(trial_num)==0)];
                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==37&filters.result(trial_num)==0)))];
            case 'positive_neutral'
                if isnan(filters.opt1_win_amount(trial_num))&&~isnan(filters.opt1_loss_amount(trial_num))
                    align_times = [align_times; single_trial_times(trial_codes==37&filters.result(trial_num)==0)];
                    trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==37&filters.result(trial_num)==0)))];
                else
                    continue
                end
            case 'negative_neutral'
                if ~isnan(filters.opt1_win_amount(trial_num))&&isnan(filters.opt1_loss_amount(trial_num))
                    align_times = [align_times; single_trial_times(trial_codes==37&filters.result(trial_num)==0)];
                    trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==37&filters.result(trial_num)==0)))];
                else
                    continue
                end
            case 'positive'
                if isnan(filters.opt1_win_amount(trial_num))&&~isnan(filters.opt1_loss_amount(trial_num))
                    align_times = [align_times; single_trial_times(trial_codes==37&filters.result(trial_num)==0)];
                    trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==37&filters.result(trial_num)==0)))];
                else
                    align_times = [align_times; single_trial_times(trial_codes==37&filters.result(trial_num)>0)];
                    trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==37&filters.result(trial_num)>0)))];
                end
            case 'negative'
                if ~isnan(filters.opt1_win_amount(trial_num))&&isnan(filters.opt1_loss_amount(trial_num))
                    align_times = [align_times; single_trial_times(trial_codes==37&filters.result(trial_num)==0)];
                    trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==37&filters.result(trial_num)==0)))];
                else
                    align_times = [align_times; single_trial_times(trial_codes==37&filters.result(trial_num)<0)];
                    trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==37&filters.result(trial_num)<0)))];
                end
            case 'baseline'
                align_times = [align_times; single_trial_times(trial_codes==1)];
                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==1)))];
            case 'left_screen'
                for fn = 1:numel(trial_codes)
                    switch trial_codes(fn)
                        case 18
                            if filters.horiz_vert(trial_num)==1
                                if filters.opt1_win_amount_pos(trial_num)<=12
                                    align_times = [align_times; single_trial_times(fn)];
                                    trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                                else
                                    continue
                                end
                            else
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            end
                        case 19
                            if filters.horiz_vert(trial_num)==1
                                if filters.opt1_win_prob_pos(trial_num)<=12
                                    align_times = [align_times; single_trial_times(fn)];
                                    trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                                else
                                    continue
                                end
                            else
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            end
                        case 20
                            if filters.horiz_vert(trial_num)==1
                                if filters.opt2_win_amount_pos(trial_num)<=22
                                    align_times = [align_times; single_trial_times(fn)];
                                    trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                                else
                                    continue
                                end
                            else
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            end
                        case 21
                            if filters.horiz_vert(trial_num)==1
                                if filters.opt2_win_prob_pos(trial_num)<=22
                                    align_times = [align_times; single_trial_times(fn)];
                                else
                                    continue
                                end
                            else
                                align_times = [align_times; single_trial_times(fn)];
                            end
                        case 22
                            if filters.horiz_vert(trial_num)==1
                                if filters.opt1_loss_amount_pos(trial_num)<=12
                                    align_times = [align_times; single_trial_times(fn)];
                                    trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                                else
                                    continue
                                end
                            else
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            end
                        case 23
                            if filters.horiz_vert(trial_num)==1
                                if filters.opt1_loss_prob_pos(trial_num)<=12
                                    align_times = [align_times; single_trial_times(fn)];
                                    trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                                else
                                    continue
                                end
                            else
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            end
                        case 24
                            if filters.horiz_vert(trial_num)==1
                                if filters.opt2_loss_amount_pos(trial_num)<=22
                                    align_times = [align_times; single_trial_times(fn)];
                                    trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                                else
                                    continue
                                end
                            else
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            end
                        case 25
                            if filters.horiz_vert(trial_num)==1
                                if filters.opt2_loss_prob_pos(trial_num)<=22
                                    align_times = [align_times; single_trial_times(fn)];
                                    trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                                else
                                    continue
                                end
                            else
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            end
                        case 26
                            if filters.horiz_vert(trial_num)==1
                                if filters.opt3_win_amount_pos(trial_num)<=32
                                    align_times = [align_times; single_trial_times(fn)];
                                    trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                                else
                                    continue
                                end
                            else
                                continue
                            end
                        case 27
                            if filters.horiz_vert(trial_num)==1
                                if filters.opt3_win_prob_pos(trial_num)<=32
                                    align_times = [align_times; single_trial_times(fn)];
                                else
                                    continue
                                end
                            else
                                continue
                            end
                        case 28
                            if filters.horiz_vert(trial_num)==1
                                if filters.opt3_loss_amount_pos(trial_num)<=32
                                    align_times = [align_times; single_trial_times(fn)];
                                else
                                    continue
                                end
                            else
                                continue
                            end
                        case 29
                            if filters.horiz_vert(trial_num)==1
                                if filters.opt3_loss_prob_pos(trial_num)<=32
                                    align_times = [align_times; single_trial_times(fn)];
                                    trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                                else
                                    continue
                                end
                            else
                                continue
                            end
                        case 30
                            if filters.horiz_vert(trial_num)==1
                                if filters.opt4_win_amount_pos(trial_num)<=42
                                    align_times = [align_times; single_trial_times(fn)];
                                    trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                                else
                                    continue
                                end
                            else
                                continue
                            end
                        case 31
                            if filters.horiz_vert(trial_num)==1
                                if filters.opt4_win_prob_pos(trial_num)<=42
                                    align_times = [align_times; single_trial_times(fn)];
                                else
                                    continue
                                end
                            else
                                continue
                            end
                        case 32
                            if filters.horiz_vert(trial_num)==1
                                if filters.opt4_loss_amount_pos(trial_num)<=42
                                    align_times = [align_times; single_trial_times(fn)];
                                    trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                                else
                                    continue
                                end
                            else
                                continue
                            end
                        case 33
                            if filters.horiz_vert(trial_num)==1
                                if filters.opt4_loss_prob_pos(trial_num)<=42
                                    align_times = [align_times; single_trial_times(fn)];
                                    trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                                else
                                    continue
                                end
                            else
                                continue
                            end
                        otherwise
                            continue
                    end
                end
            case 'right_screen'
                for fn = 1:numel(trial_codes)
                    switch trial_codes(fn)
                        case 18
                            if filters.horiz_vert(trial_num)==1
                                if filters.opt1_win_amount_pos(trial_num)>12
                                    align_times = [align_times; single_trial_times(fn)];
                                    trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                                else
                                    continue
                                end
                            else
                                continue
                            end
                        case 19
                            if filters.horiz_vert(trial_num)==1
                                if filters.opt1_win_prob_pos(trial_num)>12
                                    align_times = [align_times; single_trial_times(fn)];
                                    trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                                else
                                    continue
                                end
                            else
                                continue
                            end
                        case 20
                            if filters.horiz_vert(trial_num)==1
                                if filters.opt2_win_amount_pos(trial_num)>22
                                    align_times = [align_times; single_trial_times(fn)];
                                    trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                                else
                                    continue
                                end
                            else
                                continue
                            end
                        case 21
                            if filters.horiz_vert(trial_num)==1
                                if filters.opt2_win_prob_pos(trial_num)>22
                                    align_times = [align_times; single_trial_times(fn)];
                                    trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                                else
                                    continue
                                end
                            else
                                continue
                            end
                        case 22
                            if filters.horiz_vert(trial_num)==1
                                if filters.opt1_loss_amount_pos(trial_num)>12
                                    align_times = [align_times; single_trial_times(fn)];
                                    trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                                else
                                    continue
                                end
                            else
                                continue
                            end
                        case 23
                            if filters.horiz_vert(trial_num)==1
                                if filters.opt1_loss_prob_pos(trial_num)>12
                                    align_times = [align_times; single_trial_times(fn)];
                                    trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                                else
                                    continue
                                end
                            else
                                continue
                            end
                        case 24
                            if filters.horiz_vert(trial_num)==1
                                if filters.opt2_loss_amount_pos(trial_num)>22
                                    align_times = [align_times; single_trial_times(fn)];
                                    trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                                else
                                    continue
                                end
                            else
                                continue
                            end
                        case 25
                            if filters.horiz_vert(trial_num)==1
                                if filters.opt2_loss_prob_pos(trial_num)>22
                                    align_times = [align_times; single_trial_times(fn)];
                                    trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                                else
                                    continue
                                end
                            else
                                continue
                            end
                        case 26
                            if filters.horiz_vert(trial_num)==1
                                if filters.opt3_win_amount_pos(trial_num)>32
                                    align_times = [align_times; single_trial_times(fn)];
                                    trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                                else
                                    continue
                                end
                            else
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            end
                        case 27
                            if filters.horiz_vert(trial_num)==1
                                if filters.opt3_win_prob_pos(trial_num)>32
                                    align_times = [align_times; single_trial_times(fn)];
                                    trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                                else
                                    continue
                                end
                            else
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            end
                        case 28
                            if filters.horiz_vert(trial_num)==1
                                if filters.opt3_loss_amount_pos(trial_num)>32
                                    align_times = [align_times; single_trial_times(fn)];
                                    trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                                else
                                    continue
                                end
                            else
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            end
                        case 29
                            if filters.horiz_vert(trial_num)==1
                                if filters.opt3_loss_prob_pos(trial_num)>32
                                    align_times = [align_times; single_trial_times(fn)];
                                    trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                                else
                                    continue
                                end
                            else
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            end
                        case 30
                            if filters.horiz_vert(trial_num)==1
                                if filters.opt4_win_amount_pos(trial_num)>42
                                    align_times = [align_times; single_trial_times(fn)];
                                    trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                                else
                                    continue
                                end
                            else
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            end
                        case 31
                            if filters.horiz_vert(trial_num)==1
                                if filters.opt4_win_prob_pos(trial_num)>42
                                    align_times = [align_times; single_trial_times(fn)];
                                    trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                                else
                                    continue
                                end
                            else
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            end
                        case 32
                            if filters.horiz_vert(trial_num)==1
                                if filters.opt4_loss_amount_pos(trial_num)>42
                                    align_times = [align_times; single_trial_times(fn)];
                                    trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                                else
                                    continue
                                end
                            else
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            end
                        case 33
                            if filters.horiz_vert(trial_num)==1
                                if filters.opt4_loss_prob_pos(trial_num)>42
                                    align_times = [align_times; single_trial_times(fn)];
                                    trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                                else
                                    continue
                                end
                            else
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            end
                        otherwise
                            continue
                    end
                end
            case 'likely_outcome'
                if filters.choice(trial_num)==1
                    if filters.result(trial_num)>0
                        if filters.opt1_win_prob(trial_num)>0.5
                            align_times = [align_times; single_trial_times(trial_codes==37)];
                            trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==37)))];
                        else
                            continue
                        end
                    elseif filters.result(trial_num)<0
                        if filters.opt1_loss_prob(trial_num)>0.5
                            align_times = [align_times; single_trial_times(trial_codes==37)];
                            trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==37)))];
                        else
                            continue
                        end
                    else
                        if nansum([filters.opt1_win_prob(trial_num) filters.opt1_loss_prob(trial_num)])<0.5
                            align_times = [align_times; single_trial_times(trial_codes==37)];
                            trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==37)))];
                        else
                            continue
                        end
                    end
                elseif filters.choice(trial_num)==2
                    if filters.result(trial_num)>0
                        if filters.opt2_win_prob(trial_num)>0.5
                            align_times = [align_times; single_trial_times(trial_codes==37)];
                            trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==37)))];
                        else
                            continue
                        end
                    elseif filters.result(trial_num)<0
                        if filters.opt2_loss_prob(trial_num)>0.5
                            align_times = [align_times; single_trial_times(trial_codes==37)];
                            trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==37)))];
                        else
                            continue
                        end
                    else
                        if nansum([filters.opt2_win_prob(trial_num) filters.opt2_loss_prob(trial_num)])<0.5
                            align_times = [align_times; single_trial_times(trial_codes==37)];
                            trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==37)))];
                        else
                            continue
                        end
                    end
                elseif filters.choice(trial_num)==3
                    if filters.result(trial_num)>0
                        if filters.opt3_win_prob(trial_num)>0.5
                            align_times = [align_times; single_trial_times(trial_codes==37)];
                            trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==37)))];
                        else
                            continue
                        end
                    elseif filters.result(trial_num)<0
                        if filters.opt3_loss_prob(trial_num)>0.5
                            align_times = [align_times; single_trial_times(trial_codes==37)];
                            trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==37)))];
                        else
                            continue
                        end
                    else
                        if nansum([filters.opt3_win_prob(trial_num) filters.opt3_loss_prob(trial_num)])<0.5
                            align_times = [align_times; single_trial_times(trial_codes==37)];
                            trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==37)))];
                        else
                            continue
                        end
                    end
                elseif filters.choice(trial_num)==4
                    if filters.result(trial_num)>0
                        if filters.opt4_win_prob(trial_num)>0.5
                            align_times = [align_times; single_trial_times(trial_codes==37)];
                            trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==37)))];
                        else
                            continue
                        end
                    elseif filters.result(trial_num)<0
                        if filters.opt4_loss_prob(trial_num)>0.5
                            align_times = [align_times; single_trial_times(trial_codes==37)];
                            trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==37)))];
                        else
                            continue
                        end
                    else
                        if nansum([filters.opt4_win_prob(trial_num) filters.opt4_loss_prob(trial_num)])<0.5
                            align_times = [align_times; single_trial_times(trial_codes==37)];
                            trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==37)))];
                        else
                            continue
                        end
                    end
                else
                    continue
                end
            case 'unlikely_outcome'
                if filters.choice(trial_num)==1
                    if filters.result(trial_num)>0
                        if filters.opt1_win_prob(trial_num)<0.5
                            align_times = [align_times; single_trial_times(trial_codes==37)];
                            trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==37)))];
                        else
                            continue
                        end
                    elseif filters.result(trial_num)<0
                        if filters.opt1_loss_prob(trial_num)<0.5
                            align_times = [align_times; single_trial_times(trial_codes==37)];
                            trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==37)))];
                        else
                            continue
                        end
                    else
                        if nansum([filters.opt1_win_prob(trial_num) filters.opt1_loss_prob(trial_num)])>0.5
                            align_times = [align_times; single_trial_times(trial_codes==37)];
                            trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==37)))];
                        else
                            continue
                        end
                    end
                elseif filters.choice(trial_num)==2
                    if filters.result(trial_num)>0
                        if filters.opt2_win_prob(trial_num)<0.5
                            align_times = [align_times; single_trial_times(trial_codes==37)];
                            trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==37)))];
                        else
                            continue
                        end
                    elseif filters.result(trial_num)<0
                        if filters.opt2_loss_prob(trial_num)<0.5
                            align_times = [align_times; single_trial_times(trial_codes==37)];
                            trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==37)))];
                        else
                            continue
                        end
                    else
                        if nansum([filters.opt2_win_prob(trial_num) filters.opt2_loss_prob(trial_num)])>0.5
                            align_times = [align_times; single_trial_times(trial_codes==37)];
                            trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==37)))];
                        else
                            continue
                        end
                    end
                elseif filters.choice(trial_num)==3
                    if filters.result(trial_num)>0
                        if filters.opt3_win_prob(trial_num)<0.5
                            align_times = [align_times; single_trial_times(trial_codes==37)];
                            trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==37)))];
                        else
                            continue
                        end
                    elseif filters.result(trial_num)<0
                        if filters.opt3_loss_prob(trial_num)<0.5
                            align_times = [align_times; single_trial_times(trial_codes==37)];
                            trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==37)))];
                        else
                            continue
                        end
                    else
                        if nansum([filters.opt3_win_prob(trial_num) filters.opt3_loss_prob(trial_num)])>0.5
                            align_times = [align_times; single_trial_times(trial_codes==37)];
                            trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==37)))];
                        else
                            continue
                        end
                    end
                elseif filters.choice(trial_num)==4
                    if filters.result(trial_num)>0
                        if filters.opt4_win_prob(trial_num)<0.5
                            align_times = [align_times; single_trial_times(trial_codes==37)];
                            trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==37)))];
                        else
                            continue
                        end
                    elseif filters.result(trial_num)<0
                        if filters.opt4_loss_prob(trial_num)<0.5
                            align_times = [align_times; single_trial_times(trial_codes==37)];
                            trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==37)))];
                        else
                            continue
                        end
                    else
                        if nansum([filters.opt4_win_prob(trial_num) filters.opt4_loss_prob(trial_num)])>0.5
                            align_times = [align_times; single_trial_times(trial_codes==37)];
                            trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==37)))];
                        else
                            continue
                        end
                    end
                else
                    continue
                end
            case 'high_mag_att'
                for fn = 1:numel(trial_codes)
                    switch trial_codes(fn)
                        case 18
                            if filters.opt1_win_amount(trial_num)>0.5
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            else
                                continue
                            end
                        case 19
                            if filters.opt1_win_prob(trial_num)>0.5
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            else
                                continue
                            end
                        case 20
                            if filters.opt2_win_amount(trial_num)>0.5
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            else
                                continue
                            end
                        case 21
                            if filters.opt2_win_prob(trial_num)>0.5
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            else
                                continue
                            end
                        case 22
                            if abs(filters.opt1_loss_amount(trial_num))>0.5
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            else
                                continue
                            end
                        case 23
                            if filters.opt1_loss_prob(trial_num)>0.5
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            else
                                continue
                            end
                        case 24
                            if abs(filters.opt2_loss_amount(trial_num))>0.5
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            else
                                continue
                            end
                        case 25
                            if filters.opt2_loss_prob(trial_num)>0.5
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            else
                                continue
                            end
                        case 26
                            if filters.opt3_win_amount(trial_num)>0.5
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            else
                                continue
                            end
                        case 27
                            if filters.opt3_win_prob(trial_num)>0.5
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            else
                                continue
                            end
                        case 28
                            if abs(filters.opt3_loss_amount(trial_num))>0.5
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            else
                                continue
                            end
                        case 29
                            if filters.opt3_loss_prob(trial_num)>0.5
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            else
                                continue
                            end
                        case 30
                            if filters.opt4_win_amount(trial_num)>0.5
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            else
                                continue
                            end
                        case 31
                            if filters.opt4_win_prob(trial_num)>0.5
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            else
                                continue
                            end
                        case 32
                            if abs(filters.opt4_loss_amount(trial_num))>0.5
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            else
                                continue
                            end
                        case 33
                            if filters.opt4_loss_prob(trial_num)>0.5
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            else
                                continue
                            end
                    end
                end
            case 'low_mag_att'
                for fn = 1:numel(trial_codes)
                    switch trial_codes(fn)
                        case 18
                            if filters.opt1_win_amount(trial_num)<0.5
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            else
                                continue
                            end
                        case 19
                            if filters.opt1_win_prob(trial_num)<0.5
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            else
                                continue
                            end
                        case 20
                            if filters.opt2_win_amount(trial_num)<0.5
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            else
                                continue
                            end
                        case 21
                            if filters.opt2_win_prob(trial_num)<0.5
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            else
                                continue
                            end
                        case 22
                            if abs(filters.opt1_loss_amount(trial_num))<0.5
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            else
                                continue
                            end
                        case 23
                            if filters.opt1_loss_prob(trial_num)<0.5
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            else
                                continue
                            end
                        case 24
                            if abs(filters.opt2_loss_amount(trial_num))<0.5
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            else
                                continue
                            end
                        case 25
                            if filters.opt2_loss_prob(trial_num)<0.5
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            else
                                continue
                            end
                        case 26
                            if filters.opt3_win_amount(trial_num)<0.5
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            else
                                continue
                            end
                        case 27
                            if filters.opt3_win_prob(trial_num)<0.5
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            else
                                continue
                            end
                        case 28
                            if abs(filters.opt3_loss_amount(trial_num))<0.5
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            else
                                continue
                            end
                        case 29
                            if filters.opt3_loss_prob(trial_num)<0.5
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            else
                                continue
                            end
                        case 30
                            if filters.opt4_win_amount(trial_num)<0.5
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            else
                                continue
                            end
                        case 31
                            if filters.opt4_win_prob(trial_num)<0.5
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            else
                                continue
                            end
                        case 32
                            if abs(filters.opt4_loss_amount(trial_num))<0.5
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            else
                                continue
                            end
                        case 33
                            if filters.opt4_loss_prob(trial_num)<0.5
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            else
                                continue
                            end
                    end
                end
            case 'full_opt_info'
                for fn = 1:numel(trial_codes)
                    if fn<4
                        continue
                    end
                    switch trial_codes(fn)
                        case 18
                            if any(trial_codes(1:fn-1)==19)&&any(trial_codes(1:fn-1)==22)&&any(trial_codes(1:fn-1)==23)
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            else
                                continue
                            end
                        case 19
                            if any(trial_codes(1:fn-1)==18)&&any(trial_codes(1:fn-1)==22)&&any(trial_codes(1:fn-1)==23)
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            else
                                continue
                            end
                        case 22
                            if any(trial_codes(1:fn-1)==18)&&any(trial_codes(1:fn-1)==19)&&any(trial_codes(1:fn-1)==23)
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            else
                                continue
                            end
                        case 23
                            if any(trial_codes(1:fn-1)==18)&&any(trial_codes(1:fn-1)==19)&&any(trial_codes(1:fn-1)==22)
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            else
                                continue
                            end
                        case 20
                            if any(trial_codes(1:fn-1)==21)&&any(trial_codes(1:fn-1)==24)&&any(trial_codes(1:fn-1)==25)
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            else
                                continue
                            end
                        case 21
                            if any(trial_codes(1:fn-1)==20)&&any(trial_codes(1:fn-1)==24)&&any(trial_codes(1:fn-1)==25)
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            else
                                continue
                            end
                        case 24
                            if any(trial_codes(1:fn-1)==20)&&any(trial_codes(1:fn-1)==21)&&any(trial_codes(1:fn-1)==25)
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            else
                                continue
                            end
                        case 25
                            if any(trial_codes(1:fn-1)==20)&&any(trial_codes(1:fn-1)==21)&&any(trial_codes(1:fn-1)==24)
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            else
                                continue
                            end
                        case 26
                            if any(trial_codes(1:fn-1)==27)&&any(trial_codes(1:fn-1)==28)&&any(trial_codes(1:fn-1)==29)
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            else
                                continue
                            end
                        case 27
                            if any(trial_codes(1:fn-1)==26)&&any(trial_codes(1:fn-1)==28)&&any(trial_codes(1:fn-1)==29)
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            else
                                continue
                            end
                        case 28
                            if any(trial_codes(1:fn-1)==26)&&any(trial_codes(1:fn-1)==27)&&any(trial_codes(1:fn-1)==29)
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            else
                                continue
                            end
                        case 29
                            if any(trial_codes(1:fn-1)==26)&&any(trial_codes(1:fn-1)==27)&&any(trial_codes(1:fn-1)==28)
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            else
                                continue
                            end
                        case 30
                            if any(trial_codes(1:fn-1)==31)&&any(trial_codes(1:fn-1)==32)&&any(trial_codes(1:fn-1)==33)
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            else
                                continue
                            end
                        case 31
                            if any(trial_codes(1:fn-1)==30)&&any(trial_codes(1:fn-1)==32)&&any(trial_codes(1:fn-1)==33)
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            else
                                continue
                            end
                        case 32
                            if any(trial_codes(1:fn-1)==30)&&any(trial_codes(1:fn-1)==31)&&any(trial_codes(1:fn-1)==33)
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            else
                                continue
                            end
                        case 33
                            if any(trial_codes(1:fn-1)==30)&&any(trial_codes(1:fn-1)==31)&&any(trial_codes(1:fn-1)==32)
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            else
                                continue
                            end
                        otherwise
                            continue
                    end
                end
            case 'full_dual_prob_sweep'
                for fn = 1:numel(trial_codes)
                    if fn<8
                        continue
                    end
                    switch trial_codes(fn)
                        case 19
                            if any(trial_codes(1:fn-1)==21)&&any(trial_codes(1:fn-1)==23)&&any(trial_codes(1:fn-1)==25)&&any(trial_codes(1:fn-1)==27)&&any(trial_codes(1:fn-1)==29)&&any(trial_codes(1:fn-1)==31)&&any(trial_codes(1:fn-1)==33)
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            else
                                continue
                            end
                        case 21
                            if any(trial_codes(1:fn-1)==19)&&any(trial_codes(1:fn-1)==23)&&any(trial_codes(1:fn-1)==25)&&any(trial_codes(1:fn-1)==27)&&any(trial_codes(1:fn-1)==29)&&any(trial_codes(1:fn-1)==31)&&any(trial_codes(1:fn-1)==33)
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            else
                                continue
                            end
                        case 23
                            if any(trial_codes(1:fn-1)==19)&&any(trial_codes(1:fn-1)==21)&&any(trial_codes(1:fn-1)==25)&&any(trial_codes(1:fn-1)==27)&&any(trial_codes(1:fn-1)==29)&&any(trial_codes(1:fn-1)==31)&&any(trial_codes(1:fn-1)==33)
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            else
                                continue
                            end
                        case 25
                            if any(trial_codes(1:fn-1)==19)&&any(trial_codes(1:fn-1)==21)&&any(trial_codes(1:fn-1)==23)&&any(trial_codes(1:fn-1)==27)&&any(trial_codes(1:fn-1)==29)&&any(trial_codes(1:fn-1)==31)&&any(trial_codes(1:fn-1)==33)
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            else
                                continue
                            end
                        case 27
                            if any(trial_codes(1:fn-1)==19)&&any(trial_codes(1:fn-1)==21)&&any(trial_codes(1:fn-1)==23)&&any(trial_codes(1:fn-1)==25)&&any(trial_codes(1:fn-1)==29)&&any(trial_codes(1:fn-1)==31)&&any(trial_codes(1:fn-1)==33)
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            else
                                continue
                            end
                        case 29
                            if any(trial_codes(1:fn-1)==19)&&any(trial_codes(1:fn-1)==21)&&any(trial_codes(1:fn-1)==23)&&any(trial_codes(1:fn-1)==25)&&any(trial_codes(1:fn-1)==27)&&any(trial_codes(1:fn-1)==31)&&any(trial_codes(1:fn-1)==33)
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            else
                                continue
                            end
                        case 31
                            if any(trial_codes(1:fn-1)==19)&&any(trial_codes(1:fn-1)==21)&&any(trial_codes(1:fn-1)==23)&&any(trial_codes(1:fn-1)==25)&&any(trial_codes(1:fn-1)==27)&&any(trial_codes(1:fn-1)==29)&&any(trial_codes(1:fn-1)==33)
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            else
                                continue
                            end
                        case 33
                            if any(trial_codes(1:fn-1)==19)&&any(trial_codes(1:fn-1)==21)&&any(trial_codes(1:fn-1)==23)&&any(trial_codes(1:fn-1)==25)&&any(trial_codes(1:fn-1)==27)&&any(trial_codes(1:fn-1)==29)&&any(trial_codes(1:fn-1)==31)
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            else
                                continue
                            end
                        otherwise
                            continue
                    end
                end
            case 'full_win_prob_sweep'
                for fn = 1:numel(trial_codes)
                    if fn<4
                        continue
                    end
                    switch trial_codes(fn)
                        case 19
                            if any(trial_codes(1:fn-1)==21)&&any(trial_codes(1:fn-1)==27)&&any(trial_codes(1:fn-1)==31)
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            else
                                continue
                            end
                        case 21
                            if any(trial_codes(1:fn-1)==19)&&any(trial_codes(1:fn-1)==27)&&any(trial_codes(1:fn-1)==31)
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            else
                                continue
                            end
                        case 27
                            if any(trial_codes(1:fn-1)==19)&&any(trial_codes(1:fn-1)==21)&&any(trial_codes(1:fn-1)==31)
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            else
                                continue
                            end
                        case 31
                            if any(trial_codes(1:fn-1)==19)&&any(trial_codes(1:fn-1)==21)&&any(trial_codes(1:fn-1)==27)
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            else
                                continue
                            end
                        otherwise
                            continue
                    end
                end
            case 'full_loss_prob_sweep'
                for fn = 1:numel(trial_codes)
                    if fn<4
                        continue
                    end
                    switch trial_codes(fn)
                        case 23
                            if any(trial_codes(1:fn-1)==25)&&any(trial_codes(1:fn-1)==29)&&any(trial_codes(1:fn-1)==33)
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            else
                                continue
                            end
                        case 25
                            if any(trial_codes(1:fn-1)==23)&&any(trial_codes(1:fn-1)==29)&&any(trial_codes(1:fn-1)==33)
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            else
                                continue
                            end
                        case 29
                            if any(trial_codes(1:fn-1)==23)&&any(trial_codes(1:fn-1)==25)&&any(trial_codes(1:fn-1)==33)
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            else
                                continue
                            end
                        case 33
                            if any(trial_codes(1:fn-1)==23)&&any(trial_codes(1:fn-1)==25)&&any(trial_codes(1:fn-1)==29)
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            else
                                continue
                            end
                        otherwise
                            continue
                    end
                end
            case 'full_win_amount_sweep'
                for fn = 1:numel(trial_codes)
                    if fn<4
                        continue
                    end
                    switch trial_codes(fn)
                        case 18
                            if any(trial_codes(1:fn-1)==20)&&any(trial_codes(1:fn-1)==26)&&any(trial_codes(1:fn-1)==30)
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            else
                                continue
                            end
                        case 20
                            if any(trial_codes(1:fn-1)==18)&&any(trial_codes(1:fn-1)==26)&&any(trial_codes(1:fn-1)==30)
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            else
                                continue
                            end
                        case 26
                            if any(trial_codes(1:fn-1)==18)&&any(trial_codes(1:fn-1)==20)&&any(trial_codes(1:fn-1)==30)
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            else
                                continue
                            end
                        case 30
                            if any(trial_codes(1:fn-1)==18)&&any(trial_codes(1:fn-1)==20)&&any(trial_codes(1:fn-1)==26)
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            else
                                continue
                            end
                        otherwise
                            continue
                    end
                end
            case 'full_loss_amount_sweep'
                for fn = 1:numel(trial_codes)
                    if fn<4
                        continue
                    end
                    switch trial_codes(fn)
                        case 22
                            if any(trial_codes(1:fn-1)==24)&&any(trial_codes(1:fn-1)==28)&&any(trial_codes(1:fn-1)==32)
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            else
                                continue
                            end
                        case 24
                            if any(trial_codes(1:fn-1)==22)&&any(trial_codes(1:fn-1)==28)&&any(trial_codes(1:fn-1)==32)
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            else
                                continue
                            end
                        case 28
                            if any(trial_codes(1:fn-1)==22)&&any(trial_codes(1:fn-1)==24)&&any(trial_codes(1:fn-1)==32)
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            else
                                continue
                            end
                        case 32
                            if any(trial_codes(1:fn-1)==22)&&any(trial_codes(1:fn-1)==24)&&any(trial_codes(1:fn-1)==28)
                                align_times = [align_times; single_trial_times(fn)];
                                trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(fn)))];
                            else
                                continue
                            end
                        otherwise
                            continue
                    end
                end
            case 'top_option_tap'
                if filters.horiz_vert(trial_num)==1
                    align_times = [align_times; single_trial_times(trial_codes==2|trial_codes==3|trial_codes==6|trial_codes==7)];
                    trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==2|trial_codes==3|trial_codes==6|trial_codes==7)))];
                else
                    continue
                end
            case 'bottom_option_tap'
                if filters.horiz_vert(trial_num)==1
                    align_times = [align_times; single_trial_times(trial_codes==4|trial_codes==5|trial_codes==8|trial_codes==9)];
                    trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==4|trial_codes==5|trial_codes==8|trial_codes==9)))];
                else
                    continue
                end
            case 'left_option_tap'
                if filters.horiz_vert(trial_num)==0
                    align_times = [align_times; single_trial_times(trial_codes==2|trial_codes==3|trial_codes==6|trial_codes==7)];
                    trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==2|trial_codes==3|trial_codes==6|trial_codes==7)))];
                else
                    continue
                end
            case 'right_option_tap'
                if filters.horiz_vert(trial_num)==0
                    align_times = [align_times; single_trial_times(trial_codes==4|trial_codes==5|trial_codes==8|trial_codes==9)];
                    trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==4|trial_codes==5|trial_codes==8|trial_codes==9)))];
                else
                    continue
                end
            case 'top_left_tap'
                if ~isnan(filters.opt3_win_amount_pos(trial_num)), continue; end
                if filters.opt1_win_amount(trial_num)==11
                    align_times = [align_times; single_trial_times(trial_codes==2)];
                    trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==2)))];
                elseif filters.opt1_win_prob_pos(trial_num)==11
                    align_times = [align_times; single_trial_times(trial_codes==3)];
                    trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==3)))];
                elseif filters.opt1_loss_amount_pos(trial_num)==11
                    align_times = [align_times; single_trial_times(trial_codes==6)];
                    trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==6)))];
                elseif filters.opt1_loss_prob_pos(trial_num)==11
                    align_times = [align_times; single_trial_times(trial_codes==7)];
                    trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==7)))];
                else
                    continue
                end
            case 'top_right_tap'
                if ~isnan(filters.opt3_win_amount(trial_num)), continue; end
                if filters.horiz_vert(trial_num)==1
                    if filters.opt1_win_amount_pos(trial_num)==12
                        align_times = [align_times; single_trial_times(trial_codes==2)];
                        trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==2)))];
                    elseif filters.opt1_win_prob_pos(trial_num)==12
                        align_times = [align_times; single_trial_times(trial_codes==3)];
                        trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==3)))];
                    elseif filters.opt1_loss_amount_pos(trial_num)==12
                        align_times = [align_times; single_trial_times(trial_codes==6)];
                        trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==6)))];
                    elseif filters.opt1_loss_prob_pos(trial_num)==12
                        align_times = [align_times; single_trial_times(trial_codes==7)];
                        trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==7)))];
                    end
                elseif filters.horiz_vert(trial_num)==0
                    if filters.opt2_win_amount_pos(trial_num)==21
                        align_times = [align_times; single_trial_times(trial_codes==4)];
                        trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==4)))];
                    elseif filters.opt2_win_prob_pos(trial_num)==21
                        align_times = [align_times; single_trial_times(trial_codes==5)];
                        trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==5)))];
                    elseif filters.opt2_loss_amount_pos(trial_num)==21
                        align_times = [align_times; single_trial_times(trial_codes==8)];
                        trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==8)))];
                    elseif filters.opt2_loss_prob_pos(trial_num)==21
                        align_times = [align_times; single_trial_times(trial_codes==9)];
                        trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==9)))];
                    end
                end
            case 'bottom_left_tap'
                if ~isnan(filters.opt3_win_amount(trial_num)), continue; end
                if filters.horiz_vert(trial_num)==1
                    if filters.opt2_win_amount_pos(trial_num)==21
                        align_times = [align_times; single_trial_times(trial_codes==4)];
                        trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==4)))];
                    elseif filters.opt2_win_prob_pos(trial_num)==21
                        align_times = [align_times; single_trial_times(trial_codes==5)];
                        trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==5)))];
                    elseif filters.opt2_loss_amount_pos(trial_num)==21
                        align_times = [align_times; single_trial_times(trial_codes==8)];
                        trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==8)))];
                    elseif filters.opt2_loss_prob_pos(trial_num)==21
                        align_times = [align_times; single_trial_times(trial_codes==9)];
                        trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==9)))];
                    end
                elseif filters.horiz_vert(trial_num)==0
                    if filters.opt1_win_amount_pos(trial_num)==12
                        align_times = [align_times; single_trial_times(trial_codes==2)];
                        trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==2)))];
                    elseif filters.opt1_win_prob_pos(trial_num)==12
                        align_times = [align_times; single_trial_times(trial_codes==3)];
                        trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==3)))];
                    elseif filters.opt1_loss_amount_pos(trial_num)==12
                        align_times = [align_times; single_trial_times(trial_codes==6)];
                        trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==6)))];
                    elseif filters.opt1_loss_prob_pos(trial_num)==12
                        align_times = [align_times; single_trial_times(trial_codes==7)];
                        trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==7)))];
                    end
                end
            case 'bottom_right_tap'
                if ~isnan(filters.opt3_win_amount(trial_num)), continue; end
                if filters.opt2_win_amount_pos(trial_num)==22
                    align_times = [align_times; single_trial_times(trial_codes==4)];
                    trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==4)))];
                elseif filters.opt2_win_prob_pos(trial_num)==22
                    align_times = [align_times; single_trial_times(trial_codes==5)];
                    trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==5)))];
                elseif filters.opt2_loss_amount_pos(trial_num)==22
                    align_times = [align_times; single_trial_times(trial_codes==8)];
                    trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==8)))];
                elseif filters.opt2_loss_prob_pos(trial_num)==22
                    align_times = [align_times; single_trial_times(trial_codes==9)];
                    trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==9)))];
                end
            case 'first_inspection'
                fii = find(trial_codes>17&trial_codes<34,1);
                if ~isempty(fii)
                    align_times = [align_times; single_trial_times(fii)];
                    trial_numbers = [trial_numbers; trial_number];
                end
            case 'second_inspection'
                f2ii = find(trial_codes>17&trial_codes<34,2);
                if numel(f2ii)>=2
                    sii = f2ii(2);
                    align_times = [align_times; single_trial_times(sii)];
                    trial_numbers = [trial_numbers; trial_number];
                end
            case 'third_inspection'
                f3ii = find(trial_codes>17&trial_codes<34,3);
                if numel(f3ii)>=3
                    tii = f3ii(3);
                    align_times = [align_times; single_trial_times(tii)];
                    trial_numbers = [trial_numbers; trial_number];
                end
            case 'fourth_inspection'
                f4ii = find(trial_codes>17&trial_codes<34,4);
                if numel(f4ii)>=4
                    fii = f4ii(4);
                    align_times = [align_times; single_trial_times(fii)];
                    trial_numbers = [trial_numbers; trial_number];
                end
            otherwise
                if ischar(alignment)
                    disp('PROBLEM: FUNCTION get_align_times.m NOT PREPARED TO HANDLE THIS ALIGNMENT');
                else
                    for fn = 1:numel(trial_codes)
                        align_times = [align_times; single_trial_times(trial_codes==alignment)];
                        trial_numbers = [trial_numbers; trial_number*ones(size(single_trial_times(trial_codes==alignment)))];
                    end
                end
        end
    end