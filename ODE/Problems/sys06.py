#!/usr/bin/env python3
# -*- coding : utf-8 -*-


import numpy as np
import sys
sys.path.append('../../')

from  ODE.Solvers.Euler import euler
from  ODE.Solvers.CentDiff import centdiff
from  ODE.Rhs.rhs import Rhs
import matplotlib.pyplot as plt
from datetime import datetime
from ODE.qualityPlot.qualityplot import *
def main():

    start_time = datetime.now()
    a=4;
    c=1;
    
#y1' = (2-.5*y2)*y1
#y2' = (-1+.5*y1)*y2
#Given initial conditions
#The initial rabbit population is 6, the initial fox population is 2:

#y1(0) = 6
#y2(0) = 2
           
    # differential problem PREY-PREDATOR

    eq1 = lambda t,u : (2-0.5*u[1])*u[0]         #a*(u[0]-u[0]*u[1]);
    eq2 = lambda t,u : (-1+0.5*u[0])*u[1]         #-c*(u[1]-u[0]*u[1]);
        
    func1 = np.array([eq1,eq2])
    
    y0      = np.array([6.,2.])
    system1 = Rhs(func1, 0,15,y0,1200 )

    fwdeuler = euler.Explicit(system1)
    fet,feu  = fwdeuler.solve()   

    bwdeuler = euler.Backward(system1)
    bet,beu  = bwdeuler.solve()
    
    sieuler = euler.SemiImplicit(system1)
    siet,sieu  = sieuler.solve()
    
    Impeuler = euler.Implicit(system1)
    iet,ieu  = Impeuler.solve()
    
    midpoint = centdiff.MidPoint(system1)
    cdt,cdu  = midpoint.solve()
    
    #leapfrog    = multistep.LeapFrog(system1, 'lf1.dat')
    #lft,lfu   = leapfrog.solve()
    
#    cranknicholson = rungekutta.CrankNicholson(problem1 , 'cn1.dat')
#    cnt,cnu    = cranknicholson.solve()
#    
#    am5 = adamsmethods.AdamsMoulton(problem1, 'AM5_1.dat')
#    amt,amu = am5._5th.solve()
#        
#    rk3 = rungekutta.RK3(system1,'s1RK4.dat')
#    rkt,rku = rk3.solve()
#    
#    t0 = -10. 
#    tf = 10.      
#    u0 = y0   
#    dt = 20/100 
#    t,u = t0,u0 
#    
#    with open('BWEuler_p1.dat', 'w') as f:
#         while True:
#             f.write('%10.4f %14.10f \n' %(t, u[0]) )
#             u = euler.Backward.step(func1,t,u,dt)
#             t += dt
#             if t > tf:
#                break
#             
#    bet = np.genfromtxt('BWEuler_p1.dat',usecols=(0,)) 
#    beu = np.genfromtxt('BWEuler_p1.dat',usecols=(1,)) 
#    
#    t0 = -10. 
#    tf = 10.      
#    u0 = y0   
#    dt = 20/200 
#    t,u = t0,u0 
#    
#    with open('ImpEuler_p1.dat', 'w') as f:
#         while True:
#             f.write('%10.4f %14.10f \n' %(t, u[0]) )
#             u = euler.Implicit.step(func1,t,u,dt)
#             t += dt
#             if t > tf:
#                break
#             
#    bt = np.genfromtxt('ImpEuler_p1.dat',usecols=(0,)) 
#    bu = np.genfromtxt('ImpEuler_p1.dat',usecols=(1,)) 
#    
#    t0 = -10. 
#    tf = 10.      
#    u0 = y0   
#    dt = 20/400 
#    t,u = t0,u0 
#    with open('SIEuler_p1.dat', 'w') as f:
#         while True:
#             f.write('%10.4f %14.10f \n' %(t, u[0]) )
#             u = euler.SemiImplicit.step(func1,t,u,dt)
#             t += dt
#             if t > tf:
#                break
#             
#    siet = np.genfromtxt('SIEuler_p1.dat',usecols=(0,)) 
#    sieu = np.genfromtxt('SIEuler_p1.dat',usecols=(1,)) 
#   
#    t0 = -10. 
#    tf = 10.      
#    u0 = y0   
#    dt = 20/2000 
#    t,u = t0,u0 
#    with open('AdamsMoulton5_p1.dat', 'w') as f:
#        while True:
#            f.write('%f %f \n' %(t, u[0]) )
#            u = adamsmethods.AdamsMoulton._5th.step(func1,t,u,dt)
#            t += dt
#            if t > tf:
#               break 
#    
#    amt = np.genfromtxt('AdamsMoulton5_p1.dat',usecols=(0,)) 
#    amu = np.genfromtxt('AdamsMoulton5_p1.dat',usecols=(1,)) 









# : a*(u[0]-u[0]*u[1]);
#    eq2 = lambda t,u : -c*(u[1]-u[0]*u[1]);
#-----------------------------------------------------------------------------------------------------   
    end_time = datetime.now()
    print('Duration {}'.format(end_time - start_time))
            #**{'color': 'paraview'}
##############################################################################################    
    #fig,axs = qualityplot.Palatino(figsize=(9.5,4.5),**{'scheme':'vega'})(1,1)# ,**{'scheme':'nb'})
    
#    fig,axs = qualityplot.Elsevier(**{'scheme':'nb','grid.linestyle':'-','grid.dashes':(5,5)})(1,1,figsize=(9.5,4.5))# ,
    #fig,axs = qualityplot.Elsevier(**{'scheme':'nb','grid.linestyle':'-','grid.dashes':(5,5)})(1,1,figsize=(9.5,4.5))# 
    #fig,axs = qualityplot.StandardImproved(**{'scheme':'vega','grid.dashes':(5,5), 'linestyle':'ls1'})(1,1,figsize=(9.5,4.5))# 
#    fig,axs = qualityplot.StandardImproved(**{'scheme':'vega','grid.dashes':(5,5), 'linestyle':'paper'})(1,1,figsize=(9.5,4.5))# 
    #fig,axs = qualityplot.TexMathpazo(**{'scheme':'vega','grid.dashes':(5,5), 'linestyle':'paper'})(1,1,figsize=(9.5,4.5))# 
    #fig,axs = qualityplot.TexMathpazoBeamer(**{'scheme':'nb','grid.dashes':(5,5), 'linestyle':'paper'})(1,1,figsize=(9.5,4.5))# 
    #fig,axs = qualityplot.Beamer(**{'scheme':'nb','grid.dashes':(5,5), 'linestyle':'paper'})(1,1,figsize=(9.5,4.5))#
    #fig,axs = qualityplot.TexPaper(**{'scheme':'nb', 'linestyle':'paper'})(1,1,figsize=(9.5,4.5))# 
    #fig,axs = qualityplot.Helvet(**{'scheme':'nb', 'linestyle':'paper'})(1,1,figsize=(9.5,4.5))# 
    #fig,axs = qualityplot.HelvetBeamer(**{'linestyle':'paper'})(1,1) #,figsize=(9.5,4.5))#

    #fig,axs = qualityplot.Times(**{'linestyle':'paper'})(1,1,figsize=(12.5,4.5))# ,**{'scheme':'nb'})
    #fig,axs = qualityplot.TimesItalic(**{'linestyle':'paper'})(1,1,figsize=(9.5,4.5))#

    #fig,axs = qualityplot.TimesItalic(**{'linestyle':'paper'})(1,1,figsize=(9.5,4.5))# 
    fig,axs = PalatinoItalic()(1,1,(12.0,4.0))#
    #fig,axs = qualityplot.SansItalic()(1,1,(12.0,4.0))#
#    fig,axs = PalatinoItalicDark()(1,1)#
    #fig,axs = qualityplot.PalatinoItalicDark()(1,1)#
    
    #axs.plot(x,y,label='Analitycal', marker='o',linestyle='',markersize=4.5,color='C0',alpha=0.7)
    axs.plot(fet,feu,label='Exp Euler')
    #axs.plot(iet,ieu,label='Imp Euler')
    #axs.plot(siet,sieu,label='Sem.Imp Euler') #, linestyle=':')
    axs.plot(bet,beu,label='NR-BWDEuler') #, linestyle=':')
    #axs.plot(lft,lfu,label='Leap Frog')
    #axs.plot(cnt,cnu,label='Crank Nicholson')
    #axs.plot(amt,amu,label='Adams Moulton 5th') #,linestyle=':')
    #axs.plot(rkt,rku,label='Runge Kutta 3th')
    axs.plot(cdt,cdu,label='CentDiff 2nd')
    
    chartBox = axs.get_position()    
    legend = leg = axs.legend()
    legend.get_frame().set_linewidth(0.7)
    axs.set_position([chartBox.x0, chartBox.y0, chartBox.width*0.8, chartBox.height])
    plt.setp(legend.get_texts(), color='#2E3436')
    plt.legend(bbox_to_anchor=(1.03, 0.75), loc=2, borderaxespad=0.)
    


    #legend = leg = axs.legend()
    #legend.get_frame().set_linewidth(0.7)
    #plt.setp(legend.get_texts(), color='#2E3436')
    axs.set(xlabel = 'time' ,ylabel='y(t)')
    axs.set_title('dy$_1$/dt= -4 (y$_1$ - y$_1$y$_2$)  ;  dy$_2$/dt = -1(y$_2$ - y$_1$y$_2$)',y=1.05)
    plt.savefig('../Results/systemPrayPredator_1.pdf')   
    plt.show()


if __name__ == '__main__':
    main()









