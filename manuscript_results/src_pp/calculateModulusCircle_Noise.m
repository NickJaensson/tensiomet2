function [DataStruc] = calculateModulusCircle_Noise(DataStruc, NumSimuls, VarList)

sigma_circlefit = zeros(1,NumSimuls);
sigma_ref = zeros(1,NumSimuls);
Kmod = zeros(1,NumSimuls);

for i = 1:size(VarList,1)
    for j = 1:size(VarList, 2)
        strain = DataStruc.(VarList{i,j}).Params_phys.frac;
        strainmeasure = DataStruc.(VarList{i,j}).Params_phys.strainmeasure;

        for kk = 1:NumSimuls
            sigma_circlefit(kk) = DataStruc.(VarList{i,j}).DefSt_inv{kk}.SigmaCircleFit_dimal;
            sigma_ref(kk) = DataStruc.(VarList{i,j}).RefSt_inv{kk}.SigmaCircleFit_dimal;
            
            % Calculate dilatational modulus from prescribed strain
            if strcmp(strainmeasure,'generic')
                Kmod(kk) = (sigma_circlefit(kk) - sigma_ref(kk))/(strain - 1);
            elseif strcmp(strainmeasure,'hookean')
                Kmod(kk) = (sigma_circlefit(kk) - sigma_ref(kk))/(strain - 1);
            elseif strcmp(strainmeasure,'hencky') || strcmp(strainmeasure,'balemans')
                Kmod(kk) = (sigma_circlefit(kk) - sigma_ref(kk))/log(strain);
            elseif strcmp(strainmeasure,'pepicelli')
                Kmod(kk) = (sigma_circlefit(kk) - sigma_ref(kk))/(log(strain)/strain);
            else % Hookean
                Kmod(kk) = (sigma_circlefit(kk) - sigma_ref(kk))/(strain - 1);
            end

            DataStruc.(VarList{i,j}).DefSt_inv{kk}.KmodCircleFit_dimal = Kmod(kk);
        end

        DataStruc.(VarList{i,j}).DefSt_inv{NumSimuls+1}.KmodCircleFit_dimal_Avg = mean(Kmod);
        DataStruc.(VarList{i,j}).DefSt_inv{NumSimuls+1}.KmodCircleFit_dimal_Std = std(Kmod);

    end
end

end