
# Inputs of script
# e.g. python3 script_wholeCellModel.py rho_ratio s_ext Nrepeat Nm input_s input_p input_T input_M input_R 
# rho_ratio: occupancy of cell, i.e. ratio of cytosolic volume occupied by dry mass
# s_ext: concentration substrate s in the environment; unit: µM
# input_s: initial substrate concentration; unit: log(copy-number-per-cubic-micron)
# input_p: initial precursor concentration; unit: log(copy-number-per-cubic-micron)
# input_T: initial transporter T per cell; unit: log(copy-number-per-cell)
# input_M: initial enzyme M concentration; unit: log(copy-number-per-cubic-micron)
# input_R: initial ribosome R concentration; unit: log(copy-number-per-cubic-micron)
# Nrepeat: number of attempts to solve the problem
# Nm: number of steps in the metabolic pathway



# Outputs of script
# The script outputs a table 16 columns; the number of rows is Nrepeat, i.e. the number of attempts to solve for the optimal growth rate given the circumstances
# Outputs of python script The script outputs a table 16 columns; the number of rows is Nrepeat, i.e. the number of attempts to solve for the optimal growth rate given the circumstances The different columns correspond to: index, substrate concentration in the environment, occupancy rho, optimal growth rate found in the current attempt, K_M^* of metabolic reaction, K_M^* of ribosomal reaction, dummy output, log-concentraion** of substrate s, log-concentration** of precursor p, log-copy-number of transporter T per cell, log-concentration** of metabolic enzyme M, log-concentration** of ribosome R, volume fracton of substrate s, volume fracton of prevursor p, dummay output, volume fracton of metabolic enzyme M, volume fracton of ribosome R
# **unit of concentraion: copy number per cubic micron


import numpy as np
from scipy.optimize import minimize
import scipy.optimize
import random
import sys

rho_ratio = float(sys.argv[1]);
s_ext = float(sys.argv[2]);
Nrepeat = int(sys.argv[3]);
Nm = float(sys.argv[4]);
input_x = [];
for ii in range(5,len(sys.argv)):
    input_x.append(float(sys.argv[ii]))
input_s,input_p,input_T,input_M,input_R = input_x;

ttheta = 2.3;
ttolerance_KM = 1e-5;  # tolerance of solver
ttolerance_solver = 1e-10;  # tolerance of solver

NenergyCharge = 0;
#Nrepeat = 50;
#Nm = 10.0;
cost_p_onR = 0.0;

sint_vol0 = 1.6464e-28;  # volume of a substrate molecule s; metric unit
p_vol0 = 2.6808e-18;  # volume of a precursor molecule p; metric unit
E_vol0 = 5.7906e-26;  # volume of an enzyme E; metric unit
R_vol0 = 4.1888e-24;  # volume of a ribosome R; metric unit
tRNA_p_vol0 = 5.7906e-26;  # volume of a tertiary unit tRNA-p; metric unit
mRNA_vol0 = 10*9.2028e-27;  # this variable is not relevent to the current model

tmp_vol0 = np.array([sint_vol0,p_vol0,E_vol0,R_vol0,tRNA_p_vol0,mRNA_vol0]);
sint_vol,p_vol_dummy,E_vol,R_vol,tRNA_p_vol,mRNA_vol = tmp_vol0;

#l_area = 2e-19;
#T_area = 5e-17;

E_l = 300;  # number of monomer units in an enzyme
R_l = 7459;  # number of monomer units in a ribosome
tRNA_l = 100;  # number of monomer units in a tRNA
mRNA_l = 1000;  # this variable is not relevent to the current model

KM0_M = 130*602;  # unit: µM
KM0_R = 120*602;  # unit: µM
KM0_T = 1000;  # unit: µM
KM0_mRNA = 8*120;  # this variable is not relevent to the current model

Kcat_T = 13.7;  # unit: per second
Kcat_M = 13.7;  # unit: per second
Kcat_L = 13.7;  # unit: per second
Kcat_S = 13.7;  # unit: per second
Kcat_C = 4;  # unit: per second
Kcat_R = 22;  # unit: per second
Kcat_mRNA = 5000;  # this variable is not relevent to the current model

vol_cell = 1e-18;  # volume of a cell is 1 cubic micron, metric unit

rho_core = 0.00;  # this variable is not relevent to the current model
rho_mRNA = 0.00;  # this variable is not relevent to the current model
N_core = vol_cell*rho_core/(E_vol+sint_vol);  # this variable is not used in the current model
N_mRNA = vol_cell*rho_mRNA/mRNA_vol;  # this variable is not relevent to the current model

KM0 = np.array([KM0_M,KM0_R,KM0_mRNA]);  

tmpV_met = [sint_vol,tRNA_p_vol];
tmpV_rxn = [E_vol,R_vol,mRNA_vol];
tmpV_D = [E_vol+sint_vol,R_vol+tRNA_p_vol,mRNA_vol+R_vol];
tmpV_core = [E_vol,sint_vol];
x_V = np.array(tmpV_met+tmpV_rxn+tmpV_D+tmpV_core)  # corresponds to V_i in Eq. S3
x_H = (x_V/4.0*3.0/np.pi)**(1.0/3.0);  # corresponds to r_i in Eq. S3, i.e. radius of different molec
x_S = 4*np.pi*x_H*x_H;  # corresponds to S_i in Eq. S3

r_sint,r_p,r_M,r_R,r_mRNA,r_DM,r_DR,r_DmRNA,r_E,r_sintb = x_H.tolist();
Nmet = len(tmpV_met);
Nrxn = len(tmpV_rxn);

arr_vp = [1e3,1e4,1e5];
vol_x = np.array([sint_vol,tRNA_p_vol,0,E_vol,R_vol]); # the volume of a molecule of substrate s, precuror p, transporter T, enzyme E, ribosome R
tmp_var = 3.0;
bound_var = 6.0;
bnds = [];  # define the upper and lower bound of the concentration of substrate s, precuror p, transporter T, enzyme E, ribosome R; unit: log(copy-number-per-cubic-micron)
for ii in range(len(input_x)-1):
    ttmp = [0,0];
    if vol_x[ii]>0:
        ttmp[0] = input_x[ii]-bound_var;
        if ttmp[0]<0: ttmp[0]=0;
        #ttmp[1] = np.log(vol_cell*(rho_ratio-rho_core) / vol_x[ii]);
        ttmp[1] = input_x[ii]+bound_var;
    else:
        ttmp[0] = 0
        ttmp[1] = 50;
    bnds.append(ttmp);

# calculate the hydrodynamic radius, which is from: Trovato, Diffusion within the Cytoplasm: A Mesoscale Model of Interacting Macromolecules, 2014
def f_sp_exp2(x,r1,r2):
    #epsi = 0.51e-9*(0.25/x);
    epsi = 0.51e-9;    
    Rh = 42e-9;
    aa = 0.53;
    if r1>r2:
        rp = r2;
    else:
        rp = r1;
    rrh = (rp+1.4e-10)*1.3;
    gg = ((epsi/Rh)**2 + (epsi/rp)**2)**(-0.5*aa) / 0.22;
    output = np.exp(-1.0*gg*x);
    return output;

dict_cur_KM = {};
for jj_repeat in range(Nrepeat):
    while True:
        def f_gen_xx(input_x):
            tmp_x = np.zeros(len(input_x))
            totvol = 0;
            for ii in range(len(input_x)):
                tmplb = input_x[ii]-tmp_var;
                if tmplb<0:
                    tmplb = 0
                tmp_x[ii] = random.uniform(tmplb, input_x[ii]+tmp_var)
                totvol = totvol + np.exp(tmp_x[ii])*vol_x[ii];
            tmp_factor = np.log((rho_ratio-rho_core-rho_mRNA)*vol_cell/totvol);
            for ii in range(len(input_x)):
                if vol_x[ii]>0:
                    tmp_x[ii] = tmp_x[ii] + tmp_factor;
            iia = range(len(input_x)-1);
            out_x = tmp_x[iia];
            return out_x;
        def f_toString(xx):
            return np.array2string(xx,precision=8)
        def f_R(xx):
            model_sint,model_p,model_T,model_M = (np.exp(xx)).tolist();
            nR = ((rho_ratio-rho_core-rho_mRNA)*vol_cell - model_sint*sint_vol - (model_p)*tRNA_p_vol - (model_M)*E_vol) / R_vol
            return nR;
        def f_get_complex(xx,cur_KM):
            model_sint,model_p,model_T,model_M = (np.exp(xx)).tolist();
            nR = f_R(xx);
            cur_KM_M,cur_KM_R,cur_KM_mRNA = cur_KM.tolist();
            tmp_D_M = (model_M/Nm) / (1 + cur_KM_M/(model_sint/Nm));
            tmp_D_R = nR / (1 + cur_KM_R/(model_p));
            tmp_D_mRNA = N_mRNA / (1 + cur_KM_mRNA/nR);  
            return tmp_D_M,tmp_D_R,tmp_D_mRNA;
        def f_get_KM(xx):
            if sum(np.isnan(xx))>0: return np.exp(40)*KM0,np.exp(40)*KM0,np.exp(40)*KM0;
            txtxx = f_toString(xx);
            if txtxx in dict_cur_KM:
                return dict_cur_KM[txtxx]
            model_sint,model_p,model_T,model_M = (np.exp(xx)).tolist();
            model_R = f_R(xx)
            cur_KM = KM0.copy()
            while True:
                cur_KM_M,cur_KM_R,cur_KM_mRNA = cur_KM.tolist()
                D_M_cur,D_R_cur,D_mRNA_cur = f_get_complex(xx,cur_KM);
                D_M_cur = D_M_cur*Nm;
                x_sint_cur,x_p_cur,x_M_cur,x_R_cur = model_sint,model_p,model_M,model_R;
                
                x_M_cur = x_M_cur - D_M_cur;
                x_R_cur = x_R_cur - D_R_cur - D_mRNA_cur;
                x_mRNA_cur = N_mRNA - D_mRNA_cur;
                
                x_sint_cur = x_sint_cur - D_M_cur;
                x_p_cur = x_p_cur - D_R_cur;
                
                x_count = np.array([x_sint_cur,x_p_cur,x_M_cur,x_R_cur,x_mRNA_cur,D_M_cur,D_R_cur,D_mRNA_cur,N_core,N_core]);
                
                avg_V = np.sum(x_count*x_V)/vol_cell
                avg_S = np.sum(x_count*x_S)/vol_cell;
                avg_H = np.sum(x_count*x_H)/vol_cell;
                avg_H2 = np.sum(x_count*x_H*x_H)/vol_cell;
                avg_1 = np.sum(x_count)/vol_cell;
                log_gamma = -1.0*np.log(1-avg_V) + (x_H*avg_S + x_S*avg_H + x_V*avg_1)/(1-avg_V) + ((x_H*x_H)*(avg_S*avg_S) + 2*x_V*avg_H*avg_S)/(2*(1.0-avg_V)**2) + x_V*avg_H2*(avg_S**2)/(3*(1-avg_V)**3);
                
                lg_sint,lg_p,lg_M,lg_R,lg_mRNA,lg_DM,lg_DR,lg_DmRNA,lg_coreE,lg_coreSint = log_gamma;
                
                log_Gamma_KM_M = lg_sint + lg_M - lg_DM;
                log_Gamma_KM_R = lg_p + lg_R - lg_DR;
                log_Gamma_KM_mRNA = lg_mRNA + lg_R - lg_DmRNA;
                log_Gamma = np.array([log_Gamma_KM_M,log_Gamma_KM_R,log_Gamma_KM_mRNA]);
                
                ffactor_radius_KM_M = 1.0/r_M + 1.0/r_sint;
                ffactor_radius_KM_R = 1.0/r_p + 1.0/r_R;
                ffactor_radius_KM_mRNA = 1.0/r_mRNA + 1.0/r_R;       
                ffactor_radius = np.array([ffactor_radius_KM_M,ffactor_radius_KM_R,ffactor_radius_KM_mRNA]);
                
                pair_radius = [];
                pair_radius.append([r_M,r_sint]);
                pair_radius.append([r_p,r_R]);
                pair_radius.append([r_R,r_mRNA]);
                
                pre_KM = cur_KM.copy();
                termm_denominator = np.zeros(Nrxn);
                termm_numerator = np.zeros(Nrxn);
                for qq in range(Nrxn):
                   termm_denominator[qq] = (1+ttheta) * np.exp(log_Gamma[qq])*f_sp_exp2(avg_V,pair_radius[qq][0],pair_radius[qq][1]);
                   termm_numerator[qq] = np.exp(log_Gamma[qq])+ttheta*f_sp_exp2(avg_V,pair_radius[qq][0],pair_radius[qq][1]);
                tmp_KM = np.zeros(len(cur_KM));
                for zzz in range(len(tmp_KM)):
                    if termm_denominator[zzz]==0:
                        tmp_KM[zzz] = np.exp(40);
                    elif ((KM0[zzz] / termm_denominator[zzz] * termm_numerator[zzz])>np.exp(40)) or (np.isinf(KM0[zzz] / termm_denominator[zzz] * termm_numerator[zzz])):
                        tmp_KM[zzz] = np.exp(40);
                    else:
                        tmp_KM[zzz] = KM0[zzz] / termm_denominator[zzz] * termm_numerator[zzz];
                for qq in range(Nrxn):
                    if (np.isinf(termm_denominator[qq])):
                        tmp_exp = f_sp_exp2(avg_V,pair_radius[qq][0],pair_radius[qq][1]);
                        cur_KM[qq] = KM0[qq] / tmp_exp / (1+ttheta)
                    else:
                        cur_KM[qq] = tmp_KM[qq];
                ttmp_max = (np.abs(cur_KM-pre_KM)/pre_KM).max();
                if ttmp_max<ttolerance_KM:
                    dict_cur_KM[txtxx] = [cur_KM,log_Gamma,avg_V]
                    return cur_KM,log_Gamma,avg_V
        def f_flux_R(xx):
            if sum(np.isnan(xx))>0: return 0;
            cur_KM,log_Gamma,avg_V = f_get_KM(xx)
            cur_KM_M,cur_KM_R,cur_KM_mRNA = cur_KM.tolist()
            model_sint,model_p,model_T,model_M = (np.exp(xx)).tolist();
            nR = f_R(xx);
            rate_elongation = Kcat_R * nR / (1 + cur_KM_R/(model_p));
            return rate_elongation;
        def f_mu(xx):
            if sum(np.isnan(xx))>0: return 0;
            model_sint,model_p,model_T,model_M = (np.exp(xx)).tolist();
            model_R = f_R(xx);
            log_R = np.log(model_R);
            log_sint,log_p,log_T,log_M = xx.tolist();
            input_s,input_p,input_T,input_M,input_R
            if log_T<input_T-bound_var:
                return 0;
            elif log_R<input_R-bound_var:
                return 0;
            elif log_M<input_M-bound_var:
                return 0;
            elif log_sint<input_s-bound_var:
                return 0; 
            cur_KM,log_Gamma,avg_V = f_get_KM(xx)
            cur_KM_M,cur_KM_R,cur_KM_mRNA = cur_KM.tolist()
            return f_flux_R(xx) / ((model_M+N_core+model_T)*E_l + f_R(xx)*R_l + cost_p_onR*model_p);    
        def Obj_rule(xx):
            if np.sum(np.isnan(xx))>0: 
                return 0;
            else:
                return -1 * f_mu(xx)
        def f_flux_b(xx,indexx):
            cur_KM,log_Gamma,avg_V = f_get_KM(xx)
            cur_KM_M,cur_KM_R,cur_KM_mRNA = cur_KM.tolist()
            model_sint,model_p,model_T,model_M = (np.exp(xx)).tolist();
            output = [0,0,0,0,0,0]
            if indexx==0:
                output[indexx]  = Kcat_T*model_T/(1+KM0_T/s_ext);
            elif indexx==3:
                output[indexx] = f_flux_R(xx)
            elif indexx==5:
                output[indexx] = Kcat_M * (model_M/Nm) / (1+cur_KM_M/(model_sint/Nm));
            return output[indexx]
        
        x0 = f_gen_xx(input_x);                
        ineq_cons_lb = {'type': 'ineq', 'fun' : lambda xx: np.log(f_flux_R(xx))-(input_x[Nrxn]-bound_var)};
        ineq_cons_ub = {'type': 'ineq', 'fun' : lambda xx: (input_x[Nrxn]+bound_var)-np.log(f_flux_R(xx))};
        eq_cons_sint = {'type': 'eq', 'fun': lambda xx: (f_flux_b(xx,0)-f_flux_b(xx,5)-f_mu(xx)*np.exp(xx[0]))/np.exp(xx[0])};
        eq_cons_p = {'type': 'eq', 'fun': lambda xx: (f_flux_b(xx,5)-(NenergyCharge+1)*f_flux_R(xx)-f_mu(xx)*np.exp(xx[1]))/np.exp(xx[1])};
        cons = (ineq_cons_lb,ineq_cons_ub,eq_cons_sint,eq_cons_p);

        res = minimize(Obj_rule, x0, method='SLSQP', jac=False, constraints=cons, options={'ftol': ttolerance_solver, 'disp': True, 'maxiter' : 5000, 'disp' : False},bounds=bnds, tol=ttolerance_solver);
        
        if res.success and (f_R(res.x)>1) and ((-1*res.fun)>0):
            cur_KM,log_Gamma,avg_V = f_get_KM(res.x)
            arr_cur_KM = cur_KM.tolist()
            
            model_sint,model_p,model_T,model_M = (np.exp(res.x)).tolist();
            vector_x = np.array([model_sint,model_p,model_T,model_M,f_R(res.x)]);
            
            vector_vol = vector_x*vol_x;
            vector_logx = np.log(vector_x);
            
            arr_logx = vector_logx.tolist();
            arr_vol = vector_vol.tolist();
            headingg = [jj_repeat,s_ext,rho_ratio,-1*res.fun];
            
            arrraay = headingg;
            arrraay.extend(arr_cur_KM)
            arrraay.extend(arr_logx)
            arrraay.extend(arr_vol)
            toPrint = "\t".join(format(x, 'E') for x in arrraay);
            print(toPrint);
            sys.stdout.flush();
            break;

