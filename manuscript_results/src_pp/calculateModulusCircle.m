function [DataStruc] = calculateModulusCircle(DataStruc, VarList)

for i = 1:size(VarList,1)
    for j = 1:size(VarList, 2)

        sigma_circlefit = DataStruc.(VarList{i,j}).DefSt_inv.SigmaCircleFit_dimal;
        sigma_ref = DataStruc.(VarList{i,j}).RefSt_inv.SigmaCircleFit_dimal;
        strain = DataStruc.(VarList{i,j}).Params_phys.frac;
        strainmeasure = DataStruc.(VarList{i,j}).Params_phys.strainmeasure;

        % Calculate dilatational modulus from prescribed strain
            if strcmp(strainmeasure,'generic')
                Kmod = (sigma_circlefit - sigma_ref)/(strain - 1);
            elseif strcmp(strainmeasure,'hookean')
                Kmod = (sigma_circlefit - sigma_ref)/(strain - 1);
            elseif strcmp(strainmeasure,'hencky') || strcmp(strainmeasure,'balemans')
                Kmod = (sigma_circlefit - sigma_ref)/log(strain);
            elseif strcmp(strainmeasure,'pepicelli')
                Kmod = (sigma_circlefit - sigma_ref)/(log(strain)/strain);
            else % Hookean
                Kmod = (sigma_circlefit - sigma_ref)/(strain - 1);
            end

            DataStruc.(VarList{i,j}).DefSt_inv.KmodCircleFit_dimal = Kmod;

    end
end

end