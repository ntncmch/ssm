from Ccoder import Ccoder
import unittest
import copy
import json
import os
import shutil

class TestCcoder(unittest.TestCase):

    def setUp(self):
        dpkgRoot = os.path.join('..' ,'examples', 'noise')
        dpkg = json.load(open(os.path.join(dpkgRoot, 'ssm.json')))
        
        self.m_noise = Ccoder(dpkgRoot, copy.deepcopy(dpkg))
        
        m_diff = copy.deepcopy(dpkg)    
        del m_diff['reactions'][2]['white_noise']
        del m_diff['reactions'][3]['white_noise']
        m_diff['sde'] = {
            'drift': [
                {'name': 'r0_nyc', 'f': 0.0, 'transformation': 'log(r0_nyc)'},
                {'name': 'r0_paris', 'f': 0.0, 'transformation': 'log(r0_paris)'}
            ],
            'dispertion': [['vol',0],[0,'vol']]
        }
        self.m_diff = Ccoder(dpkgRoot, m_diff)

        m_noise2 = copy.deepcopy(dpkg)    
        del m_noise2['reactions'][2]['white_noise']
        del m_noise2['reactions'][3]['white_noise']
        m_noise2['reactions'][0]['white_noise'] = {'name':'noise_SI', 'sd': 'sto'}
        m_noise2['reactions'][1]['white_noise'] = {'name':'noise_SI2', 'sd': 'sto'}
        self.m_noise2 = Ccoder(dpkgRoot, m_noise2)

        m_noise3 = copy.deepcopy(dpkg)    
        del m_noise3['reactions'][2]['white_noise']
        del m_noise3['reactions'][3]['white_noise']
        m_noise3['reactions'][4]['white_noise'] = {'name':'noise_SI', 'sd': 'sto'}
        m_noise3['reactions'][5]['white_noise'] = {'name':'noise_SI2', 'sd': 'sto'}
        self.m_noise3 = Ccoder(dpkgRoot, m_noise3)

        m_noise4 = copy.deepcopy(dpkg)    
        del m_noise4['reactions'][2]['white_noise']
        del m_noise4['reactions'][3]['white_noise']
        m_noise4['reactions'][8]['white_noise'] = {'name':'noise_SI', 'sd': 'sto'}
        m_noise4['reactions'][9]['white_noise'] = {'name':'noise_SI2', 'sd': 'sto'}
        self.m_noise4 = Ccoder(dpkgRoot, m_noise4)

        m_noise5 = copy.deepcopy(dpkg)    
        del m_noise5['reactions'][2]['white_noise']
        del m_noise5['reactions'][3]['white_noise']
        m_noise5['reactions'][10]['white_noise'] = {'name':'noise_SI', 'sd': 'sto'}
        m_noise5['reactions'][11]['white_noise'] = {'name':'noise_SI2', 'sd': 'sto'}
        self.m_noise5 = Ccoder(dpkgRoot, m_noise5)

        m_noise6 = copy.deepcopy(dpkg)    
        m_noise6['reactions'][4]['white_noise'] = {'name':'noise_SI', 'sd': 'sto'}
        m_noise6['reactions'][5]['white_noise'] = {'name':'noise_SI2', 'sd': 'sto'}
        self.m_noise6 = Ccoder(dpkgRoot, m_noise6)

        m_noise7 = copy.deepcopy(dpkg)    
        m_noise7['reactions'][4]['white_noise'] = {'name':'noise_SI23', 'sd': 'sto'}
        m_noise7['reactions'][5]['white_noise'] = {'name':'noise_SI24', 'sd': 'sto'}
        self.m_noise7 = Ccoder(dpkgRoot, m_noise7)

        m_diff2 = copy.deepcopy(dpkg)    
        del m_diff2['reactions'][2]['white_noise']
        del m_diff2['reactions'][3]['white_noise']
        m_diff2['sde'] = copy.deepcopy(m_diff['sde'])
        m_diff2['reactions'].append({'from': 'R_paris', 'to': 'I_paris', 'rate': 'correct_rate(v)', 'description':'testing'})
        m_diff2['reactions'].append({'from': 'R_nyc', 'to': 'I_nyc', 'rate': 'correct_rate(v)', 'description':'testing'})
        self.m_diff2 = Ccoder(dpkgRoot, m_diff2)
    
        # # Tutorial with drift, to test handling of full sde
        # dpkgRoot = os.path.join('..' ,'examples', 'tutorial-with-drift')
        # dpkg = json.load(open(os.path.join(dpkgRoot, 'ssm.json')))
        # m_sde = copy.deepcopy(dpkg)
        # m_sde['sde']['drift'] = [{ 'name': 'r0', 'drift': 'r0' }]
        # m_sde['sde']['dispersion'] = [['vol*r0']]
        # self.m_sde = Ccoder(dpkgRoot, m_sde)

    def test_eval_Q(self):
        calc_Q = self.m_noise.eval_Q()

        # testing env sto only
 #       # order: I_nyc, I_paris, S_nyc, S_paris, inc_nyc, inc_all_out
        # order: S_nyc, S_paris, I_nyc, I_paris, inc_nyc, inc_all_out
        term_paris = '((((r0_paris/N_paris*v*I_paris)*S_paris)*((sto)**2))*((r0_paris/N_paris*v*I_paris)*S_paris)))'
        term_nyc = '((((r0_nyc/N_nyc*v*I_nyc)*S_nyc)*((sto)**2))*((r0_nyc/N_nyc*v*I_nyc)*S_nyc)))'
        Q_cm = calc_Q["no_dem_sto"]["Q_cm"]
        self.assertEqual(Q_cm[0][0], '((-1)*'+term_nyc+'*(-1)')
        self.assertEqual(Q_cm[0][1], 0)
        self.assertEqual(Q_cm[0][2], '((-1)*'+term_nyc+'*(1)')
        self.assertEqual(Q_cm[0][3], 0)
        self.assertEqual(Q_cm[0][4], '((-1)*'+term_nyc+'*(1)')
        self.assertEqual(Q_cm[0][5], 0)
        self.assertEqual(Q_cm[1][0], 0)
        self.assertEqual(Q_cm[1][1], '((-1)*'+term_paris+'*(-1)')
        self.assertEqual(Q_cm[1][2], 0)
        self.assertEqual(Q_cm[1][3], '((-1)*'+term_paris+'*(1)')
        self.assertEqual(Q_cm[1][4], 0)
        self.assertEqual(Q_cm[1][5], 0)
        self.assertEqual(Q_cm[2][0], '((1)*'+term_nyc+'*(-1)')
        self.assertEqual(Q_cm[2][1], 0)
        self.assertEqual(Q_cm[2][2], '((1)*'+term_nyc+'*(1)')
        self.assertEqual(Q_cm[2][3], 0)
        self.assertEqual(Q_cm[2][4], '((1)*'+term_nyc+'*(1)')
        self.assertEqual(Q_cm[2][5], 0)
        self.assertEqual(Q_cm[3][0], 0)
        self.assertEqual(Q_cm[3][1], '((1)*'+term_paris+'*(-1)')
        self.assertEqual(Q_cm[3][2], 0)
        self.assertEqual(Q_cm[3][3], '((1)*'+term_paris+'*(1)')
        self.assertEqual(Q_cm[3][4], 0)
        self.assertEqual(Q_cm[3][5], 0)
        self.assertEqual(Q_cm[4][0], '((1)*'+term_nyc+'*(-1)')
        self.assertEqual(Q_cm[4][1], 0)
        self.assertEqual(Q_cm[4][2], '((1)*'+term_nyc+'*(1)')
        self.assertEqual(Q_cm[4][3], 0)
        self.assertEqual(Q_cm[4][4], '((1)*'+term_nyc+'*(1)')
        self.assertEqual(Q_cm[4][5], 0)
        self.assertEqual(Q_cm[5][2], 0)
        self.assertEqual(Q_cm[5][3], 0)
        self.assertEqual(Q_cm[5][0], 0)
        self.assertEqual(Q_cm[5][1], 0)
        self.assertEqual(Q_cm[5][4], 0)
        self.assertEqual(Q_cm[5][5], 0)

        # testing dem sto only
        term1_p = '(mu_b_paris*N_paris))'
        term2_p = '((r0_paris/N_paris*v*I_paris)*S_paris))'
        term3_p = '((correct_rate(v))*I_paris))'
        term4_p = '((mu_d_paris)*S_paris))'
        term5_p = '((mu_d_paris)*I_paris))'
        term1_n = '(mu_b_nyc*N_nyc))'
        term2_n = '((r0_nyc/N_nyc*v*I_nyc)*S_nyc))'
        term3_n = '((correct_rate(v))*I_nyc))'
        term4_n = '((mu_d_nyc)*S_nyc))'
        term5_n = '((mu_d_nyc)*I_nyc))'
        # order: I_nyc, I_paris, S_nyc, S_paris, inc_all , inc_nyc
        Q_cm = calc_Q["no_env_sto"]["Q_cm"]
        self.assertEqual(Q_cm[0][0], '((1)*'+term1_n+'*(1) + ((-1)*'+term2_n+'*(-1) + ((-1)*'+term4_n+'*(-1)')
        self.assertEqual(Q_cm[0][1], 0)
        self.assertEqual(Q_cm[0][2], '((-1)*'+term2_n+'*(1)')
        self.assertEqual(Q_cm[0][3], 0)
        self.assertEqual(Q_cm[0][4], '((-1)*'+term2_n+'*(1)')
        self.assertEqual(Q_cm[0][5], 0)
        self.assertEqual(Q_cm[1][0], 0)
        self.assertEqual(Q_cm[1][1], '((1)*'+term1_p+'*(1) + ((-1)*'+term2_p+'*(-1) + ((-1)*'+term4_p+'*(-1)')
        self.assertEqual(Q_cm[1][2], 0)
        self.assertEqual(Q_cm[1][3], '((-1)*'+term2_p+'*(1)')
        self.assertEqual(Q_cm[1][4], 0)
        self.assertEqual(Q_cm[1][5], 0)
        self.assertEqual(Q_cm[2][0], '((1)*'+term2_n+'*(-1)')
        self.assertEqual(Q_cm[2][1], 0)
        self.assertEqual(Q_cm[2][2], '((1)*'+term2_n+'*(1) + ((-1)*'+term3_n+'*(-1) + ((-1)*'+term5_n+'*(-1)')
        self.assertEqual(Q_cm[2][3], 0)
        self.assertEqual(Q_cm[2][4], '((1)*'+term2_n+'*(1)')
        self.assertEqual(Q_cm[2][5], '((-1)*'+term3_n+'*(1) + ((-1)*'+term5_n+'*(1)')
        self.assertEqual(Q_cm[3][0], 0)
        self.assertEqual(Q_cm[3][1], '((1)*'+term2_p+'*(-1)')
        self.assertEqual(Q_cm[3][2], 0)
        self.assertEqual(Q_cm[3][3], '((1)*'+term2_p+'*(1) + ((-1)*'+term3_p+'*(-1) + ((-1)*'+term5_p+'*(-1)')
        self.assertEqual(Q_cm[3][4], 0)
        self.assertEqual(Q_cm[3][5], '((-1)*'+term3_p+'*(1) + ((-1)*'+term5_p+'*(1)')
        self.assertEqual(Q_cm[4][0], '((1)*'+term2_n+'*(-1)')
        self.assertEqual(Q_cm[4][1], 0)
        self.assertEqual(Q_cm[4][2], '((1)*'+term2_n+'*(1)')
        self.assertEqual(Q_cm[4][3], 0)
        self.assertEqual(Q_cm[4][4], '((1)*'+term2_n+'*(1)')
        self.assertEqual(Q_cm[4][5], 0)
        self.assertEqual(Q_cm[5][0], 0)
        self.assertEqual(Q_cm[5][1], 0)
        self.assertEqual(Q_cm[5][2], '((1)*'+term3_n+'*(-1) + ((1)*'+term5_n+'*(-1)')
        self.assertEqual(Q_cm[5][3], '((1)*'+term3_p+'*(-1) + ((1)*'+term5_p+'*(-1)')
        self.assertEqual(Q_cm[5][4], 0)
        self.assertEqual(Q_cm[5][5], '((1)*'+term3_n+'*(1) + ((1)*'+term3_p+'*(1) + ((1)*'+term5_n+'*(1) + ((1)*'+term5_p+'*(1)')
        
        

    # def test_eval_Q_tricky_cases(self):

    #     calc_Q = self.m_noise2.eval_Q()
    #     # testing env sto only for m_noise2 : WN on U->S
    #     term_p = '(((mu_b_paris*N_paris)*((sto)**2))*(mu_b_paris*N_paris)))'
    #     term_n = '(((mu_b_nyc*N_nyc)*((sto)**2))*(mu_b_nyc*N_nyc)))'
    #     for i in range(5):
    #         for j in range(5):
    #             if i==2 and j == 2:
    #                 self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][i][j],'((1)*'+term_n+'*(1)')
    #             elif i==3 and j == 3:
    #                 self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][i][j],'((1)*'+term_p+'*(1)')
    #             else:
    #                 self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][i][j],0)


    #     calc_Q = self.m_noise3.eval_Q()
    #     # testing env sto only for m_noise3 : WN on I->R
    #     term_p = '((((correct_rate(v))*I_paris)*((sto)**2))*((correct_rate(v))*I_paris)))'
    #     term_n = '((((correct_rate(v))*I_nyc)*((sto)**2))*((correct_rate(v))*I_nyc)))'
    #     Q_cm = calc_Q["no_dem_sto"]["Q_cm"]
    #     for i in range(5):
    #         for j in range(5):
    #             if i==0 and j == 0:
    #                 self.assertEqual(Q_cm[i][j], '((-1)*'+term_n+'*(-1)')
    #             elif i==1 and j == 1:
    #                 self.assertEqual(Q_cm[i][j], '((-1)*'+term_p+'*(-1)')
    #             elif i==0  and j == 4:
    #                 self.assertEqual(Q_cm[i][j], '((-1)*'+term_n+'*(1)')
    #             elif i==4  and j == 0:
    #                 self.assertEqual(Q_cm[i][j], '((1)*'+term_n+'*(-1)')
    #             elif i==1  and j == 4:
    #                 self.assertEqual(Q_cm[i][j], '((-1)*'+term_p+'*(1)')
    #             elif i==4  and j == 1:
    #                 self.assertEqual(Q_cm[i][j], '((1)*'+term_p+'*(-1)')
    #             elif i==4  and j == 4:
    #                 self.assertEqual(Q_cm[i][j], '((1)*'+term_p+'*(1) + ((1)*'+term_n+'*(1)')
    #             else:
    #                 self.assertEqual(Q_cm[i][j], 0)
        

        # calc_Q = self.m_noise4.eval_Q()
        # # testing env sto only for m_noise4 : WN on I->U
        # term_p = '((((mu_d_paris)*I_paris)*((sto)**2))*((mu_d_paris)*I_paris)))'
        # term_n = '((((mu_d_nyc)*I_nyc)*((sto)**2))*((mu_d_nyc)*I_nyc)))'
        # Q_cm = calc_Q["no_dem_sto"]["Q_cm"]
        # for i in range(5):
        #     for j in range(5):
        #         if i==0 and j == 0:
        #             self.assertEqual(Q_cm[i][j], '((-1)*'+term_n+'*(-1)')
        #         elif i==1 and j == 1:
        #             self.assertEqual(Q_cm[i][j], '((-1)*'+term_p+'*(-1)')
        #         elif i==1  and j == 4:
        #             self.assertEqual(Q_cm[i][j], '((-1)*'+term_p+'*(1)')
        #         elif i==4  and j == 1:
        #             self.assertEqual(Q_cm[i][j], '((1)*'+term_p+'*(-1)')
        #         elif i==0  and j == 4:
        #             self.assertEqual(Q_cm[i][j], '((-1)*'+term_n+'*(1)')
        #         elif i==4  and j == 0:
        #             self.assertEqual(Q_cm[i][j], '((1)*'+term_n+'*(-1)')
        #         elif i==4  and j == 4:
        #             self.assertEqual(Q_cm[i][j], '((1)*'+term_p+'*(1) + ((1)*'+term_n+'*(1)')
        #         else:
        #             self.assertEqual(Q_cm[i][j], 0)

                    

        # calc_Q = self.m_noise5.eval_Q()
        # # testing env sto only for m_noise5 : WN on R->U
        # for i in range(5):
        #     for j in range(5):
        #         self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][i][j],0)

        # calc_Q = self.m_noise6.eval_Q()
        # # testing env sto only for m_noise6 : correlated WN on I->R and S->I
        # term1_p = '((r0_paris/N_paris*v*I_paris)*S_paris)'
        # term2_p = '((correct_rate(v))*I_paris)'
        # term1_n = '((r0_nyc/N_nyc*v*I_nyc)*S_nyc)'
        # term2_n = '((correct_rate(v))*I_nyc)'
        # Q_cm = calc_Q["no_dem_sto"]["Q_cm"]
        # self.assertEqual(Q_cm[0][0], '((1)*(('+term1_n+'*((sto)**2))*'+term1_n+') + (-1)*(('+term2_n+'*((sto)**2))*'+term1_n+'))*(1) + ((1)*(('+term1_n+'*((sto)**2))*'+term2_n+') + (-1)*(('+term2_n+'*((sto)**2))*'+term2_n+'))*(-1)')
        # self.assertEqual(Q_cm[0][1], 0)
        # self.assertEqual(Q_cm[0][2], '((1)*(('+term1_n+'*((sto)**2))*'+term1_n+') + (-1)*(('+term2_n+'*((sto)**2))*'+term1_n+'))*(-1)')
        # self.assertEqual(Q_cm[0][3], 0)
        # self.assertEqual(Q_cm[0][4], '((1)*(('+term1_n+'*((sto)**2))*'+term2_n+') + (-1)*(('+term2_n+'*((sto)**2))*'+term2_n+'))*(1)')
        # self.assertEqual(Q_cm[0][5], '((1)*(('+term1_n+'*((sto)**2))*'+term1_n+') + (-1)*(('+term2_n+'*((sto)**2))*'+term1_n+'))*(1)')
        # self.assertEqual(Q_cm[1][0], 0)
        # self.assertEqual(Q_cm[1][1], '((1)*(('+term1_p+'*((sto)**2))*'+term1_p+') + (-1)*(('+term2_p+'*((sto)**2))*'+term1_p+'))*(1) + ((1)*(('+term1_p+'*((sto)**2))*'+term2_p+') + (-1)*(('+term2_p+'*((sto)**2))*'+term2_p+'))*(-1)')
        # self.assertEqual(Q_cm[1][2], 0)
        # self.assertEqual(Q_cm[1][3], '((1)*(('+term1_p+'*((sto)**2))*'+term1_p+') + (-1)*(('+term2_p+'*((sto)**2))*'+term1_p+'))*(-1)')
        # self.assertEqual(Q_cm[1][4], '((1)*(('+term1_p+'*((sto)**2))*'+term2_p+') + (-1)*(('+term2_p+'*((sto)**2))*'+term2_p+'))*(1)')
        # self.assertEqual(Q_cm[1][5], 0)
        # self.assertEqual(Q_cm[2][0], '((-1)*(('+term1_n+'*((sto)**2))*'+term1_n+'))*(1) + ((-1)*(('+term1_n+'*((sto)**2))*'+term2_n+'))*(-1)')
        # self.assertEqual(Q_cm[2][1], 0)
        # self.assertEqual(Q_cm[2][2], '((-1)*(('+term1_n+'*((sto)**2))*'+term1_n+'))*(-1)')
        # self.assertEqual(Q_cm[2][3], 0)
        # self.assertEqual(Q_cm[2][4], '((-1)*(('+term1_n+'*((sto)**2))*'+term2_n+'))*(1)')
        # self.assertEqual(Q_cm[2][5], '((-1)*(('+term1_n+'*((sto)**2))*'+term1_n+'))*(1)')
        # self.assertEqual(Q_cm[3][0], 0)
        # self.assertEqual(Q_cm[3][1], '((-1)*(('+term1_p+'*((sto)**2))*'+term1_p+'))*(1) + ((-1)*(('+term1_p+'*((sto)**2))*'+term2_p+'))*(-1)')
        # self.assertEqual(Q_cm[3][2], 0)
        # self.assertEqual(Q_cm[3][3], '((-1)*(('+term1_p+'*((sto)**2))*'+term1_p+'))*(-1)')
        
        # self.assertEqual(Q_cm[3][4], '((-1)*(('+term1_p+'*((sto)**2))*'+term2_p+'))*(1)')
        # self.assertEqual(Q_cm[3][5], 0)
        # self.assertEqual(Q_cm[4][0], '((1)*(('+term2_n+'*((sto)**2))*'+term1_n+'))*(1) + ((1)*(('+term2_n+'*((sto)**2))*'+term2_n+'))*(-1)')
        # self.assertEqual(Q_cm[4][1], '((1)*(('+term2_p+'*((sto)**2))*'+term1_p+'))*(1) + ((1)*(('+term2_p+'*((sto)**2))*'+term2_p+'))*(-1)')
        # self.assertEqual(Q_cm[4][2], '((1)*(('+term2_n+'*((sto)**2))*'+term1_n+'))*(-1)')
        # self.assertEqual(Q_cm[4][3], '((1)*(('+term2_p+'*((sto)**2))*'+term1_p+'))*(-1)')
        # self.assertEqual(Q_cm[4][4], '((1)*(('+term2_p+'*((sto)**2))*'+term2_p+'))*(1) + ((1)*(('+term2_n+'*((sto)**2))*'+term2_n+'))*(1)')
        # self.assertEqual(Q_cm[4][5], '((1)*(('+term2_n+'*((sto)**2))*'+term1_n+'))*(1)')
        # self.assertEqual(Q_cm[5][0], '((1)*(('+term1_n+'*((sto)**2))*'+term1_n+'))*(1) + ((1)*(('+term1_n+'*((sto)**2))*'+term2_n+'))*(-1)')
        # self.assertEqual(Q_cm[5][1], 0)
        # self.assertEqual(Q_cm[5][2], '((1)*(('+term1_n+'*((sto)**2))*'+term1_n+'))*(-1)')
        # self.assertEqual(Q_cm[5][3], 0)
        # self.assertEqual(Q_cm[5][4], '((1)*(('+term1_n+'*((sto)**2))*'+term2_n+'))*(1)')
        # self.assertEqual(Q_cm[5][5], '((1)*(('+term1_n+'*((sto)**2))*'+term1_n+'))*(1)')
        
        # calc_Q = self.m_noise7.eval_Q()
        # # testing env sto only for m_noise7 : uncorrelated WN on I->R and S->I
        # term1_p = '((r0_paris/N_paris*v*I_paris)*S_paris)'
        # term2_p = '((correct_rate(v))*I_paris)'
        # term1_n = '((r0_nyc/N_nyc*v*I_nyc)*S_nyc)'
        # term2_n = '((correct_rate(v))*I_nyc)'
        # Q_cm = calc_Q["no_dem_sto"]["Q_cm"]
        # self.assertEqual(Q_cm[0][0], '((1)*(('+term1_n+'*((sto)**2))*'+term1_n+'))*(1) + ((-1)*(('+term2_n+'*((sto)**2))*'+term2_n+'))*(-1)')
        # self.assertEqual(Q_cm[0][1], 0)
        # self.assertEqual(Q_cm[0][2], '((1)*(('+term1_n+'*((sto)**2))*'+term1_n+'))*(-1)')
        # self.assertEqual(Q_cm[0][3], 0)
        # self.assertEqual(Q_cm[0][4], '((-1)*(('+term2_n+'*((sto)**2))*'+term2_n+'))*(1)')
        # self.assertEqual(Q_cm[0][5], '((1)*(('+term1_n+'*((sto)**2))*'+term1_n+'))*(1)')
        # self.assertEqual(Q_cm[1][0], 0)
        # self.assertEqual(Q_cm[1][1], '((1)*(('+term1_p+'*((sto)**2))*'+term1_p+'))*(1) + ((-1)*(('+term2_p+'*((sto)**2))*'+term2_p+'))*(-1)')
        # self.assertEqual(Q_cm[1][2], 0)
        # self.assertEqual(Q_cm[1][3], '((1)*(('+term1_p+'*((sto)**2))*'+term1_p+'))*(-1)')
        # self.assertEqual(Q_cm[1][4], '((-1)*(('+term2_p+'*((sto)**2))*'+term2_p+'))*(1)')
        # self.assertEqual(Q_cm[1][5], 0)
        # self.assertEqual(Q_cm[2][0], '((-1)*(('+term1_n+'*((sto)**2))*'+term1_n+'))*(1)')
        # self.assertEqual(Q_cm[2][1], 0)
        # self.assertEqual(Q_cm[2][2], '((-1)*(('+term1_n+'*((sto)**2))*'+term1_n+'))*(-1)')
        # self.assertEqual(Q_cm[2][3], 0)
        # self.assertEqual(Q_cm[2][4], 0)
        # self.assertEqual(Q_cm[2][5], '((-1)*(('+term1_n+'*((sto)**2))*'+term1_n+'))*(1)')
        # self.assertEqual(Q_cm[3][0], 0)
        # self.assertEqual(Q_cm[3][1], '((-1)*(('+term1_p+'*((sto)**2))*'+term1_p+'))*(1)')
        # self.assertEqual(Q_cm[3][2], 0)
        # self.assertEqual(Q_cm[3][3], '((-1)*(('+term1_p+'*((sto)**2))*'+term1_p+'))*(-1)')
        
        # self.assertEqual(Q_cm[3][4], 0)
        # self.assertEqual(Q_cm[3][5], 0)
        # self.assertEqual(Q_cm[4][0], '((1)*(('+term2_n+'*((sto)**2))*'+term2_n+'))*(-1)')
        # self.assertEqual(Q_cm[4][1], '((1)*(('+term2_p+'*((sto)**2))*'+term2_p+'))*(-1)')
        # self.assertEqual(Q_cm[4][2], 0)
        # self.assertEqual(Q_cm[4][3], 0)
        # self.assertEqual(Q_cm[4][4], '((1)*(('+term2_p+'*((sto)**2))*'+term2_p+'))*(1) + ((1)*(('+term2_n+'*((sto)**2))*'+term2_n+'))*(1)')
        # self.assertEqual(Q_cm[4][5], 0)
        # self.assertEqual(Q_cm[5][0], '((1)*(('+term1_n+'*((sto)**2))*'+term1_n+'))*(1)')
        # self.assertEqual(Q_cm[5][1], 0)
        # self.assertEqual(Q_cm[5][2], '((1)*(('+term1_n+'*((sto)**2))*'+term1_n+'))*(-1)')
        # self.assertEqual(Q_cm[5][3], 0)
        # self.assertEqual(Q_cm[5][4], 0)
        # self.assertEqual(Q_cm[5][5], '((1)*(('+term1_n+'*((sto)**2))*'+term1_n+'))*(1)')

    def test_jac(self):
        step_ode_sde = self.m_noise.step_ode_sde()
        jac = self.m_diff.jac(step_ode_sde['sf'])

        '''
        testing Jacobians
          - Equations are processed store in the cache in the same order they are in the 'reactions' section of the ssm.json file
          - State variables are indexed in the order of appearance in reactions' definitions ("from" processed before "to"):
                { 
                    S_nyc: 0, S_paris: 1, I_nyc: 2, I_paris: 3 
                }
          - jac_obs terms correspond to the incoming rates for each observed incidence variables, for each accumulator, ordered 
            following their appearance in the reactions definitions. It is an array of cache indexes of derivatives with respect 
            to each state variable. It may contain repetitions of several accumulators feed the same state variable. Here are the 
            accumulators:
                {
                    "S_nyc to I_nyc": { 
                        accumulator_index: 0, 
                        incidence variable: {
                            'name': nyc_inc
                            'index': 0
                        }
                    },
                    "I_nyc to R_nyc": { 
                        accumulator_index: 1, 
                        incidence variable: {
                            'name': all_inc_out
                            'index': 1
                        }
                    },
                    "I_paris to R_paris": { 
                        accumulator_index: 2, 
                        incidence variable: {
                            'name': all_inc_out
                            'index': 1
                        }
                    },
                    "I_nyc to U": { 
                        accumulator_index: 3, 
                        incidence variable: {
                            'name': all_inc_out
                            'index': 2
                        }
                    },
                    "I_paris to U": { 
                        accumulator_index: 4, 
                        incidence variable: {
                            'name': all_inc_out
                            'index': 2
                        }
                    }
                }
            - jac_obs_diff is the same as jac_obs, but for diffusing variables. It is an array of cache indexes of derivatives with respect 
            to each diffusing variable. It may contain repetitions of several accumulators feed the same state variable. Here are the 
            accumulators.
        '''

        # S ode - ((r0/N*v*I)*S) - ((mu_d)*S) + (mu_b*N)
        # NYC
        self.assertEqual(jac['caches'][jac['jac'][0][0]], '-X[ORDER_I_nyc]*diffed[ORDER_diff__r0_nyc]*gsl_vector_get(par,ORDER_v)/gsl_spline_eval(calc->spline[ORDER_N_nyc],t,calc->acc[ORDER_N_nyc])-gsl_spline_eval(calc->spline[ORDER_mu_d_nyc],t,calc->acc[ORDER_mu_d_nyc])')
        self.assertEqual(jac['caches'][jac['jac'][0][1]], '0')
        self.assertEqual(jac['caches'][jac['jac'][0][2]], '-X[ORDER_S_nyc]*diffed[ORDER_diff__r0_nyc]*gsl_vector_get(par,ORDER_v)/gsl_spline_eval(calc->spline[ORDER_N_nyc],t,calc->acc[ORDER_N_nyc])')
        self.assertEqual(jac['caches'][jac['jac'][0][3]], '0')

        # Paris
        self.assertEqual(jac['caches'][jac['jac'][1][0]], '0')
        self.assertEqual(jac['caches'][jac['jac'][1][1]], '-X[ORDER_I_paris]*diffed[ORDER_diff__r0_paris]*gsl_vector_get(par,ORDER_v)/gsl_spline_eval(calc->spline[ORDER_N_paris],t,calc->acc[ORDER_N_paris])-gsl_spline_eval(calc->spline[ORDER_mu_d_paris],t,calc->acc[ORDER_mu_d_paris])')
        self.assertEqual(jac['caches'][jac['jac'][1][2]], '0')
        self.assertEqual(jac['caches'][jac['jac'][1][3]], '-X[ORDER_S_paris]*diffed[ORDER_diff__r0_paris]*gsl_vector_get(par,ORDER_v)/gsl_spline_eval(calc->spline[ORDER_N_paris],t,calc->acc[ORDER_N_paris])')

        # I ode - ((v)*I) - ((mu_d)*I) + ((r0/N*v*I)*S)
        # NYC
        self.assertEqual(jac['caches'][jac['jac'][2][0]], 'X[ORDER_I_nyc]*diffed[ORDER_diff__r0_nyc]*gsl_vector_get(par,ORDER_v)/gsl_spline_eval(calc->spline[ORDER_N_nyc],t,calc->acc[ORDER_N_nyc])')
        self.assertEqual(jac['caches'][jac['jac'][2][1]], '0')
        self.assertEqual(jac['caches'][jac['jac'][2][2]], '-gsl_spline_eval(calc->spline[ORDER_mu_d_nyc],t,calc->acc[ORDER_mu_d_nyc])-(gsl_vector_get(par,ORDER_v))+X[ORDER_S_nyc]*diffed[ORDER_diff__r0_nyc]*gsl_vector_get(par,ORDER_v)/gsl_spline_eval(calc->spline[ORDER_N_nyc],t,calc->acc[ORDER_N_nyc])')
        self.assertEqual(jac['caches'][jac['jac'][2][3]], '0')

        # Paris
        self.assertEqual(jac['caches'][jac['jac'][3][0]], '0')
        self.assertEqual(jac['caches'][jac['jac'][3][1]], 'X[ORDER_I_paris]*diffed[ORDER_diff__r0_paris]*gsl_vector_get(par,ORDER_v)/gsl_spline_eval(calc->spline[ORDER_N_paris],t,calc->acc[ORDER_N_paris])')
        self.assertEqual(jac['caches'][jac['jac'][3][2]], '0')
        self.assertEqual(jac['caches'][jac['jac'][3][3]], '-gsl_spline_eval(calc->spline[ORDER_mu_d_paris],t,calc->acc[ORDER_mu_d_paris])-(gsl_vector_get(par,ORDER_v))+X[ORDER_S_paris]*diffed[ORDER_diff__r0_paris]*gsl_vector_get(par,ORDER_v)/gsl_spline_eval(calc->spline[ORDER_N_paris],t,calc->acc[ORDER_N_paris])')        
        
        # testing jac_obs
        # nyc_inc terms (S_nyc to I_nyc accumulator): (r0/N*v*I)*S
        self.assertEqual(jac['caches'][jac['jac_obs'][0][0]], 'X[ORDER_I_nyc]*diffed[ORDER_diff__r0_nyc]*gsl_vector_get(par,ORDER_v)/gsl_spline_eval(calc->spline[ORDER_N_nyc],t,calc->acc[ORDER_N_nyc])')
        self.assertEqual(jac['caches'][jac['jac_obs'][0][1]], '0')
        self.assertEqual(jac['caches'][jac['jac_obs'][0][2]], 'X[ORDER_S_nyc]*diffed[ORDER_diff__r0_nyc]*gsl_vector_get(par,ORDER_v)/gsl_spline_eval(calc->spline[ORDER_N_nyc],t,calc->acc[ORDER_N_nyc])')
        self.assertEqual(jac['caches'][jac['jac_obs'][0][3]], '0')
        self.assertEqual(jac['caches'][jac['jac_obs_diff'][0][0]['value']], 'X[ORDER_I_nyc]*X[ORDER_S_nyc]*gsl_vector_get(par,ORDER_v)/gsl_spline_eval(calc->spline[ORDER_N_nyc],t,calc->acc[ORDER_N_nyc])')
        self.assertEqual(jac['caches'][jac['jac_obs_diff'][0][1]['value']], '0')

        # all_inc_out terms from I_nyc to R_nyc: correct_rate(v)(I_nyc + I_paris) + mu_d_nyc * I_nyc + mu_d_paris * I_paris
        self.assertEqual(jac['caches'][jac['jac_obs'][1][0]], '0')
        self.assertEqual(jac['caches'][jac['jac_obs'][1][1]], '0')
        self.assertEqual(jac['caches'][jac['jac_obs'][1][2]], 'gsl_spline_eval(calc->spline[ORDER_mu_d_nyc],t,calc->acc[ORDER_mu_d_nyc])+(gsl_vector_get(par,ORDER_v))')
        self.assertEqual(jac['caches'][jac['jac_obs'][1][3]], 'gsl_spline_eval(calc->spline[ORDER_mu_d_paris],t,calc->acc[ORDER_mu_d_paris])+(gsl_vector_get(par,ORDER_v))')
        self.assertEqual(jac['caches'][jac['jac_obs_diff'][1][0]['value']], '0')
        self.assertEqual(jac['caches'][jac['jac_obs_diff'][1][1]['value']], '0')

        # all_inc_out terms from I_paris to R_paris: correct_rate(v)(I_nyc + I_paris) + mu_d_nyc * I_nyc + mu_d_paris * I_paris
        # -> same as above
        self.assertEqual(jac['caches'][jac['jac_obs'][2][0]], '0')
        self.assertEqual(jac['caches'][jac['jac_obs'][2][1]], '0')
        self.assertEqual(jac['caches'][jac['jac_obs'][2][2]], 'gsl_spline_eval(calc->spline[ORDER_mu_d_nyc],t,calc->acc[ORDER_mu_d_nyc])+(gsl_vector_get(par,ORDER_v))')
        self.assertEqual(jac['caches'][jac['jac_obs'][2][3]], 'gsl_spline_eval(calc->spline[ORDER_mu_d_paris],t,calc->acc[ORDER_mu_d_paris])+(gsl_vector_get(par,ORDER_v))')
        self.assertEqual(jac['caches'][jac['jac_obs_diff'][2][0]['value']], '0')
        self.assertEqual(jac['caches'][jac['jac_obs_diff'][2][1]['value']], '0')

        # all_inc_out terms from I_nyc to U: correct_rate(v)(I_nyc + I_paris) + mu_d_nyc * I_nyc + mu_d_paris * I_paris
        # -> same as above
        self.assertEqual(jac['caches'][jac['jac_obs'][3][0]], '0')
        self.assertEqual(jac['caches'][jac['jac_obs'][3][1]], '0')
        self.assertEqual(jac['caches'][jac['jac_obs'][3][2]], 'gsl_spline_eval(calc->spline[ORDER_mu_d_nyc],t,calc->acc[ORDER_mu_d_nyc])+(gsl_vector_get(par,ORDER_v))')
        self.assertEqual(jac['caches'][jac['jac_obs'][3][3]], 'gsl_spline_eval(calc->spline[ORDER_mu_d_paris],t,calc->acc[ORDER_mu_d_paris])+(gsl_vector_get(par,ORDER_v))')
        self.assertEqual(jac['caches'][jac['jac_obs_diff'][3][0]['value']], '0')
        self.assertEqual(jac['caches'][jac['jac_obs_diff'][3][1]['value']], '0')

        # all_inc_out terms from I_paris to U: correct_rate(v)(I_nyc + I_paris) + mu_d_nyc * I_nyc + mu_d_paris * I_paris
        # -> same as above
        self.assertEqual(jac['caches'][jac['jac_obs'][4][0]], '0')
        self.assertEqual(jac['caches'][jac['jac_obs'][4][1]], '0')
        self.assertEqual(jac['caches'][jac['jac_obs'][4][2]], 'gsl_spline_eval(calc->spline[ORDER_mu_d_nyc],t,calc->acc[ORDER_mu_d_nyc])+(gsl_vector_get(par,ORDER_v))')
        self.assertEqual(jac['caches'][jac['jac_obs'][4][3]], 'gsl_spline_eval(calc->spline[ORDER_mu_d_paris],t,calc->acc[ORDER_mu_d_paris])+(gsl_vector_get(par,ORDER_v))')
        self.assertEqual(jac['caches'][jac['jac_obs_diff'][4][0]['value']], '0')
        self.assertEqual(jac['caches'][jac['jac_obs_diff'][4][1]['value']], '0')
        

    # def test_cache_special_function_C(self):

    #     caches = map(lambda x: self.m_diff.make_C_term(x, False), ['sin(2*PI*(t +r0))', 'sin(2*PI*(t +r0))', 'sin(2*PI*(t +r0)) + correct_rate(v)'])
    #     sf = self.m_diff.cache_special_function_C(caches)

    #     self.assertEqual(sf, ['sin(2*PI*(r0+t))', 'ssm_correct_rate(gsl_vector_get(par,ORDER_v),dt)'])
    #     self.assertEqual(caches, ['_sf[0]', '_sf[0]', '_sf[1]+_sf[0]'])


    # def test_cache_special_function_C_nested(self):

    #     caches = map(lambda x: self.m_diff.make_C_term(x, False), ['pow(correct_rate(correct_rate(v)), pow(2,4))'])
    #     sf = self.m_diff.cache_special_function_C(caches)

    #     self.assertEqual(sf, ['pow(ssm_correct_rate(ssm_correct_rate(gsl_vector_get(par,ORDER_v),dt),dt),pow(2,4))'])
    #     self.assertEqual(caches, ['_sf[0]'])


if __name__ == '__main__':
    unittest.main()
