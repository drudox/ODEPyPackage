#!/usr/bin/env python3
# -*- coding : utf-8 -*-


import numpy as np
#import rhs
#import euler
#import rungekutta
#import multistep
import matplotlib.pyplot as plt
#from datetime import datetime
import qualityplot
#mport centdiff



#    start_time = datetime.now()

def plot_sin():
    '''Plot several sinusoids'''
    n = 6
    #ax = plt.subplot(111)
    L = 2*np.pi
    x = np.linspace(0, L)
    shift = np.linspace(0, L, n, endpoint=False)
    fig,axs = qualityplot.PalatinoItalic(**{'scheme':'nb', 'lines.linewidth':'4'})(1,1)
    for s in shift:
        axs.plot(x, np.sin(x + s))
    axs.margins(0)
    axs.set_title('Sinusoids', y=1.05)
    axs.set_xlabel('x')
    axs.set_ylabel('y')
    return axs
        
           
plot_sin()

plt.show()




# differential problem 
#-----------------------------------------------------------------------------------------------------   
    #print('Duration {}'.format(end_time - start_time))
            #**{'color': 'paraview'}
    #fig,axs = qualityplot.Standard(**{'scheme':'vega'})(1,1,figsize=(8.5,4.5))# ,**{'scheme':'nb'})
    #fig,axs = qualityplot.StandardImproved(**{'scheme':'vega'} )(1,1,figsize=(9.5,4.5))# ,**{'scheme':'nb'})
    #fig,axs = qualityplot.Elsevier()(1,1,figsize=(9.5,4.5))# ,**{'scheme':'nb'})
    #fig,axs = qualityplot.Elsevier(figsize=(9.5,4.5),**{'scheme':'mycolor'})(1,1)# ,**{'scheme':'nb'})
    #fig,axs = qualityplot.TexMathpazo(figsize=(9.5,4.5))(1,1)# ,**{'scheme':'nb'})
    #fig,axs = qualityplot.TexMathpazo(figsize=(9.5,4.5),**{'scheme':'mycolor'})(1,1)# ,**{'scheme':'nb'})
    #fig,axs = qualityplot.TexMathpazoBeamer(figsize=(6.5,4.5),**{'scheme':'mycolor'})(1,1)# ,**{'scheme':'nb'})
    #fig,axs = qualityplot.TexPaper()(1,1,figsize=(8.5,4.5))# ,**{'scheme':'nb'})
    #fig,axs = qualityplot.Helvet()(1,1,figsize=(9.5,4.5))# ,**{'scheme':'nb'})
    #fig,axs = qualityplot.HelvetBeamer()(1,1,figsize=(8.5,5.5))# ,**{'scheme':'nb'})
    #fig,axs = qualityplot.Beamer()(1,1,figsize=(8.5,5.5))# ,**{'scheme':'nb'})
    #fig,axs = qualityplot.Times()(1,1,figsize=(9.5,4.5))# ,**{'scheme':'nb'})
    #fig,axs = qualityplot.TimesItalicDefault(figsize=(9.5,4.5))(1,1)# ,**{'scheme':'nb'})
    #fig,axs = qualityplot.TimesItalicImproved()(1,1,figsize=(9.5,4.5))# ,**{'scheme':'nb'})
    #fig,axs = qualityplot.TimesItalicModified()(1,1,figsize=(9.5,4.5))# ,**{'scheme':'nb'})
    #fig,axs = qualityplot.SansModified(figsize=(9.5,4.5),**{'scheme':'vega'} )(1,1)
    #fig,axs = qualityplot.Fonts(figsize=(9.5,4.5))(1,1)# ,**{'scheme':'nb'})
    #fig,axs = qualityplot.Palatino(figsize=(9.5,4.5))(1,1)# ,**{'scheme':'nb'})
    #fig,axs = qualityplot.PalatinoItalic(**{'scheme':'nb'})(1,1)# ,**{'scheme':'nb'})
    
    #axs.plot(x,y,label='Analitycal', marker='o',linestyle='',markersize=4.5,color='C0',alpha=0.7)
  #  axs.plot(fet,feu,label='Exp Euler',color='C0')
  #  #axs.plot(iet,ieu,label='Imp Euler')
  #  #axs.plot(siet,sieu,label='Sem.Imp Euler') #, linestyle=':')
  #  
  #  axs.plot(bet,beu,label='NR-BWDEuler',color='C1') #, linestyle=':')
  #  #axs.plot(rkt,rku,label='RK3')
  #  #axs.plot(lft,lfu,label='Leap Frog') #, linestyle=':')
  #  #axs.plot(lft,lfu,label='Leap Frog')
  #  #axs.plot(cnt,cnu,label='Crank Nicholson')
  #  #axs.plot(amt,amu,label='Adams Moulton 5th') #,linestyle=':')
  #  #axs.plot(rkt,rku,label='Runge Kutta 3th')
  #  #axs.plot(cdt,cdu,label='Mid Point 2nd')
  #  axs.plot(cd2t,cd2u,label='CentDiff 2nd',linestyle='-',color='C2')

  #  
  #  chartBox = axs.get_position()    
  #  legend = leg = axs.legend()
  #  legend.get_frame().set_linewidth(0.7)
  #  axs.set_position([chartBox.x0, chartBox.y0, chartBox.width*0.8, chartBox.height])
  #  plt.setp(legend.get_texts(), color='#2E3436')
  #  plt.legend(bbox_to_anchor=(1.03, 0.75), loc=2, borderaxespad=0.)
  #  
  #  axs.set(xlabel = 'time' ,ylabel='y(t)')
  #  axs.set_title('dy$_1$/dt= -4 (y$_1$ - y$_1$y$_2$)  ;  dy$_2$/dt = -1(y$_2$ - y$_1$y$_2$)',y=1.05)
  #  
  #  #plt.savefig('system02.pdf')   
  #  plt.show()


#if __name__ == '__main__':
#    main()









