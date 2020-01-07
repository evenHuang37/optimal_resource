#!/usr/bin/env python3
# -*- coding:utf-8 -*-
__author__ = 'Even Huang'
import random
import numpy as np
import matplotlib.pyplot as plt
import pylab as pl
from first.generate_ER import ER_network as ER
# from BA_wjc import BA_network as BA
import re

class SIS_model:
    def __init__(self ,N ,edge_p, beta, miu, t, p_0 ,network_type, method, repeat_times):
        self.beta = beta
        self.miu = miu
        self.p_0 = p_0
        self.t = t
        self.network_type = network_type
        self.method = method
        self.repeat_times = repeat_times
        self.N = N
        self.edge_p = edge_p


        # 选择不同的网络进行生产SIS模型
        if self.network_type == "ER_network":
            er_network = ER(self.N, self.edge_p,title = "ER network")
            self.net_mat = er_network.Create_ER_network()
            self.SIS_list = np.ones(len(self.net_mat))  # 记录SIS状态的列表
        # elif self.network=="BA_network":
        #     ba_network=BA(N=3,p=0.006,N_end=1000,m0=3,title="BA network")
        #     self.net_mat=ba_network.Create_BA_network(\
        #         ba_network.Create_ER_network())
        else:
            print("Please choose correct network")

    #获得感染密度
    def get_density(self):
         return np.sum(self.SIS_list==2)/float(len(self.net_mat))

    #初始传播
    def Generate_SIS_model(self):
        t=self.t
        net_mat=self.net_mat
        self.SIS_list = np.ones(len(self.net_mat))  # 新的一轮传播
        if self.method=="random_node":
            # 随机选择一个节点作为传播者
            rand_picked_I=np.random.choice(len(net_mat))
            self.SIS_list[rand_picked_I]=2
        elif self.method=="max_node":
            # 选择度最大的节点作为传播者
            row_sum=net_mat.sum(axis=1)
            row_sum_max=np.where(row_sum==np.max(row_sum))
            self.SIS_list[row_sum_max]=2
        # elif  re.match("density:", self.method) != None :
        elif self.method == "density":
            # 初始感染密度为p_0
            infected_num = int(self.p_0 * len(net_mat))
            for i in range(infected_num):
                infected_node = random.randint(0, len(net_mat) - 1)
                while(self.SIS_list[infected_node] == 2):
                    infected_node = random.randint(0, len(net_mat) - 1)
                self.SIS_list[infected_node] = 2
        else:
            print("Please choose a method")


        #分布记录SIS状态的列表
        SIS_t_seq_s = [np.sum(self.SIS_list==1)]
        SIS_t_seq_i = [np.sum(self.SIS_list==2)]


        for times in range(t):
            s_index = np.where(self.SIS_list == 1)#s_index为状态为s的序列
            # 在本次循环中需要根据已存在的SIS状态(old_SIS_list)进行操作，
            # 即old_SIS_list为（t-1）时刻的SIS状态
            old_SIS_list=self.SIS_list
            # I进行感染，需要对每一个s状态的节点进行遍历
            for s in s_index[0]:
                # s为状态为s的节点
                s_row=net_mat[s]# i是状态为s的一行
                link_index=np.where(s_row==1)#在i这行中存在连线的序号
                link_node=old_SIS_list[link_index]#和s这个节点相连的点的序列
                S_to_I_sum=np.sum(link_node==2)#s这个节点和状态i连线的个数
                P_of_infection=1-(1-self.beta)**S_to_I_sum#该节点被感染的概率
                if random.random()<P_of_infection:
                    self.SIS_list[s]=2
            # I可能康复
            old_I_index=np.where(old_SIS_list==2)
            for i in old_I_index[0]:
                if random.random()<self.miu:
                    self.SIS_list[i]=1
            # 记录每一时刻s,i状态的变化情况
            SIS_t_seq_s.append(np.sum(self.SIS_list==1))
            SIS_t_seq_i.append(np.sum(self.SIS_list==2))
            print(self.get_density())
        # 返回t时间内s，i的数量变化序列
        return SIS_t_seq_s,SIS_t_seq_i

    def Graw_SIS(self):
        #重复实验
        infect = np.zeros(self.t+1)
        suspect = np.zeros(self.t+1)
        for i in range(self.repeat_times):
            SIS_t_seq_s,SIS_t_seq_i= self.Generate_SIS_model()
            infect += SIS_t_seq_i
            suspect += SIS_t_seq_s
        infect = infect/self.repeat_times
        suspect = suspect/self.repeat_times
        pl.subplot(111)
        pl.plot(suspect, '-g', label='Susceptibles')
        pl.plot(infect, '-r', label='Infectious')
        # pl.plot(SIS_t_seq_r, '-k', label='Recovereds')
        pl.legend(loc=0)
        pl.title('SIS_Model based on ' + self.network_type)
        pl.xlabel('Time')
        pl.ylabel('Infectious Susceptibles and Recovereds')
        pl.xlabel('Time')
        pl.show()


    def main(self):
        self.Graw_SIS()