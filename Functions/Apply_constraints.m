function [fr_flag,P_flag,P_fr_flag,None_flag]=Apply_constraints

     fr_flag=0;          % Pick resonance frequency constraint detection mode                
     P_flag=0;           % Pick resonance Interpulse constraint detection mode             
     P_fr_flag=0;        % Pick resonance Interpulse & resonance frequency constraints detection mode
     None_flag=0;

     answer_apply_constraint = questdlg('Apply constraint?', ...
    'Apply constraint?', ...
    'Yes','No','Other');
    switch answer_apply_constraint
        case 'Yes'
            answer_Choose_constraint = questdlg('Choose constraint', ...
            'Choose constraint', ...
            'P','fr','P & fr','Other');
            % Handle response
            switch answer_Choose_constraint
            case 'P'           
                 P_flag=1;    
            case 'fr'
                 fr_flag=1;
            case 'P & fr'
                 P_fr_flag=1;
            end
        case 'No'
            None_flag=1;
    end


end