#-*- coding : utf-8 -*- 
#
#

import sys
import os
import matplotlib
import matplotlib.pyplot as plt 
import numpy as np
from matplotlib.ticker import MultipleLocator, FormatStrFormatter 
from matplotlib.axes import Axes
from cycler import cycler
from matplotlib import cm
from collections import OrderedDict
from matplotlib.ticker import AutoMinorLocator
from . basequalityplot import BasePlot

#######################################################################
#######################################################################

class Standard(BasePlot):

    def __init__(self,**kwargs):
        
        self.parameters = kwargs
        
        if 'box' not in self.parameters.keys():
            self.parameters['box'] = '#AAAAAA' 
        if 'axeslabel' not in self.parameters.keys():
            self.parameters['axeslabel'] = '#AAAAAA'
        if 'axes.linewidth' not in self.parameters.keys():
            self.parameters['axes.linewidth'] = 0.7
        if 'xtickcolor' not in self.parameters.keys():
            self.parameters['xtickcolor'] ='#AAAAAA' 
        if 'ytickcolor' not in self.parameters.keys(): 
            self.parameters['ytickcolor'] = 'gray' 
        if 'gridcolor' not in self.parameters.keys(): 
            self.parameters['gridcolor'] = 'gray' #'#AAAAAA' #'#dddddd'
        if 'font' not in self.parameters.keys(): 
            self.parameters['font'] = 'serif'
        if 'fontstyle' not in self.parameters.keys():
            self.parameters['fontstyle'] = 'italic'
        if 'fontsize' not in self.parameters.keys():
            self.parameters['fontsize'] = 10.0
        if 'legendfontsize' not in self.parameters.keys():
            self.parameters['legendfontsize'] =10.0
        if 'legendEdgeColor' not in self.parameters.keys():
            self.parameters['legendEdgeColor'] = '#AAAAAA' #'#dddddd'
        if 'scheme' not in self.parameters.keys():
            self.parameters['scheme'] = 'nb'
        if 'cycle' not in self.parameters.keys():
            self.parameters['cycle']  = self.cycle('0')
        if 'axes.linewidth' not in self.parameters.keys():
            self.parameters['axes.linewidth'] = 0.7
        if 'grid.dashes' not in self.parameters.keys():
            self.parameters['grid.dashes'] = (5,5)
        if 'grid.linestyle' not in self.parameters.keys():
            self.parameters['grid.linestyle'] = '--' 
        if 'linestyle' not in self.parameters.keys():
            self.parameters['linestyle'] = self.linestyles('paper')
        if 'cycle' not in self.parameters.keys():
            self.parameters['cycle'] = self.cycle('0')
        if 'grid.alpha' not in self.parameters.keys():
            self.parameters['grid.alpha'] = '1'
        if 'grid.linewidth' not in self.parameters.keys():
            self.parameters['grid.linewidth'] = 0.7


        self.parameters.update(kwargs)

        super().__init__(**self.parameters)

#############################################################################################
#--------------------------------------------------------------------------------------------
#       LEGEND OUTSIDE THE BOX put this lines in the main script (where the plot is called)
#--------------------------------------------------------------------------------------------
#   chartBox = axs.get_position()    
#   legend = leg = axs.legend()
#   legend.get_frame().set_linewidth(0.7)
#   axs.set_position([chartBox.x0, chartBox.y0, chartBox.width*0.8, chartBox.height])
#   plt.setp(legend.get_texts(), color='#2E3436')
#   plt.legend(bbox_to_anchor=(1.03, 0.75), loc=2, borderaxespad=0.)
#   axs.set(xlabel = 'time' ,ylabel='y(t)')
#   axs.set_title('dy$_1$/dt=-10(t-1)y$_1$',y=1.05)
#---------------------------------------------------------------------------------------            
class StandardImproved(Standard):
 
    def __init__(self,**kwargs):
        
        self.parameters = kwargs
    
        kwarg = {'label' : ['time [s]', 'y(t)= f(t)'], 'box':'#555555', 'axeslabel':'#555555', 'fontsize':'16', 'legendfontsize':'12' , 'labelsize':[16,14,14] , 'subplot.bottom':0.175 ,**kwargs}   
      
        self.parameters.update(kwarg)
       
        if 'axeslabel' not in self.parameters.keys():
            self.parameters['axeslabel'] = '#2E3436'
        if 'text.color' not in self.parameters.keys():
            self.parameters['text.color'] = 'k'
        if 'xtickcolor' not in self.parameters.keys():
            self.parameters['xtickcolor'] = 'gray' 
        if 'ytickcolor' not in self.parameters.keys():
            self.parameters['ytickcolor'] = 'gray' 
        if 'linestyle' not in self.parameters.keys():
            self.parameters['linestyle'] = self.linestyles('paper')
        if 'scheme' not in self.parameters.keys():
            self.parameters['scheme'] = 'nb' 
        if 'fontsize' not in self.parameters.keys():
            self.parameters['fontsize'] = 14.0
        if 'legendfontsize' not in self.parameters.keys():
            self.parameters['legendfontsize'] =12.0
        if 'xtick.labelsize' not in self.parameters.keys():
            self.parameters['xtick.labelsize'] = 16
        if 'ytick.labelsize' not in self.parameters.keys():
            self.parameters['ytick.labelsize'] = 16
        if 'legend.loc' not in self.parameters.keys():
            self.parameters['legend.loc'] = 'upper right'
        if 'subplot.bottom' not in self.parameters.keys():
            self.parameters['subplot.bottom'] = 0.175
        if 'subplot.right' not in self.parameters.keys():
            self.parameters['subplot.right'] = 0.93
  
        super().__init__(**self.parameters )
        
    
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
########################################################################
########################################################################

class Elsevier(StandardImproved):

    def __init__(self, **kwargs):
        
        self.parameters = kwargs
    
        #kwarg = {'label' : ['$time [s]$', '$y(t)$'], 'box':'#555555', 'axeslabel':'#555555', 'fontsize':'10', 'legendfontsize':'10' , 'labelsize':[12,12,12] ,  **kwargs}   
        
        kwarg = {'label' : ['$time [s]$', '$y(t)$'], 'box':'#555555', 'axeslabel':'#555555' ,   'subplot.bottom':0.175 ,**kwargs}   
        self.parameters.update(kwarg)
        
        self.parameters['font'] = 'Gulliver-Regular'
        if 'grid.alpha' not in self.parameters.keys():
            self.parameters['grid.alpha'] = '1'
        if 'box' not in self.parameters.keys():
            self.parameters['box'] = 'gray'
        if 'axes.linewidth' not in self.parameters.keys():
            self.parameters['axes.linewidth'] = 0.7
        if 'axeslabel' not in self.parameters.keys():
            self.parameters['axeslabel'] = '#2E3436'
        if 'text.color' not in self.parameters.keys():
            self.parameters['text.color'] = 'k'
        if 'xtickcolor' not in self.parameters.keys():
            self.parameters['xtickcolor'] = 'gray' 
        if 'ytickcolor' not in self.parameters.keys():
            self.parameters['ytickcolor'] = 'gray' 
        #if 'linestyle' not in self.parameters.keys():
        #    self.parameters['linestyle'] = self.linestyles('paper')
        if 'scheme' not in self.parameters.keys(): 
            self.parameters['scheme'] = 'vega' 
        

        super().__init__(**self.parameters )
        
#---------------------------------------------------------------------------------------            
             
class TexMathpazo(StandardImproved):

    def __init__(self, **kwargs):
        
        self.parameters = kwargs
            
        kwarg = {'package' : [r'\usepackage{amsmath}' , r'\usepackage{mathpazo}'], 'label' : ['$time [s]$', '$y(t)$', '$y\'= -10(t-1)x$'], 'text':True , 'labelsize':[16,14,14], **kwargs }   
        
        self.parameters.update(kwarg)
        
        self.parameters['fontstyle'] = 'italic'
        #self.parameters['grid.linestyle'] = '--' 
        #self.parameters['grid.alpha'] = '1'
        self.parameters['box'] = 'gray' 
        self.parameters['axes.linewidth'] = 1
        #self.parameters['axeslabel'] = 'gray'
        #self.parameters['xtickcolor'] = 'gray' 
        #self.parameters['ytickcolor'] = 'gray' 
        #self.parameters['linestyle'] = 'paper'
        if 'scheme' not in self.parameters.keys():
            self.parameters['scheme'] = 'grayscale'
        #self.parameters['fontsize'] = 12.0
        #self.parameters['legendfontsize'] =12.0
        #self.parameters['xtick.labelsize'] = 12
        #self.parameters['ytick.labelsize'] = 12
        #self.parameters['gridcolor'] = '#dddddd'
        #self.parameters['linestyle'] = 'paper'
        #self.parameters['gridcolor'] = 'gray'
        #self.parameters['grid.dashes'] = (5,5)
        #self.parameters['grid.linewidth'] = 0.7
        self.parameters['cycle']  = self.cycle('0')
        self.parameters['label'] = ['$time [s]$', '$y(t)$']
        self.parameters['legend.loc'] = 'upper right'
      
        super().__init__(**self.parameters )
 
##---------------------------------------------------------------------------------------            
##---------------------------------------------------------------------------------------            
             
class TexMathpazoBeamer(TexMathpazo):

    def __init__(self, **kwargs):
        
        self.parameters = kwargs
      
            
        kwarg = {'package' : [r'\usepackage{amsmath}' , r'\usepackage{mathpazo}'] , 'text':True , 'labelsize':[16,14,14],**kwargs }   
        self.parameters.update(kwarg)
        
        #self.parameters['box'] = 'gray' 
        #self.parameters['axes.linewidth'] = 0.7
        #self.parameters['axeslabel'] = '#2E3436'
        #self.parameters['xtickcolor'] = '#2E3436' 
        #self.parameters['ytickcolor'] = '#2E3436' 
        if 'scheme' not in self.parameters.keys():
            self.parameters['scheme'] = 'vega'
        self.parameters['fontsize'] = 16.0
        self.parameters['legendfontsize'] =14.0
        self.parameters['xtick.labelsize'] = 16
        self.parameters['ytick.labelsize'] = 16
        #self.parameters['gridcolor'] = 'gray'
        self.parameters['label'] = ['$time [s]$', '$y(t)$']
        #self.parameters['legend.loc'] = 'upper left'
        self.parameters['title.colors'] = 'k'
        #self.parameters['grid.dashes'] = (4,4)
        #self.parameters['grid.linewidth'] = 0.7
        self.parameters['subplot.bottom'] = 0.175
     #   self.parameters['subplot.right'] = 0.93
        self.parameters['lines.linewidth'] = 4
        if 'subplot.right' not in self.parameters.keys():
            self.parameters['subplot.right'] = 0.90
  
    
        super().__init__(**self.parameters )
 #----------------------------------------------------------------------------------------------------------

class Beamer(StandardImproved):
   
     def __init__(self , **kwargs):
        
        self.parameters = kwargs
         
        
        
        self.parameters['font'] = 'Neo Euler'       
        self.parameters['box'] = '#2E3436'    #'gray' 
        self.parameters['axeslabel'] = '#2E3436'
        self.parameters['fontsize'] = 18.0
        self.parameters['legendfontsize'] =14.0
        self.parameters['xtick.labelsize'] = 16
        self.parameters['ytick.labelsize'] = 16
        if 'scheme' not in self.parameters.keys():
            self.parameters['scheme'] = 'vega'
        self.parameters['legend.loc'] = 'upper right'
        self.parameters['lines.linewidth'] = 4
        self.parameters['grid.linewidth'] = 0.8
        self.parameters['title.colors'] = 'k'
        self.parameters['figuresize'] = [5,4]
        #super().get({'title.colors':'#2E3436'})
        self.parameters['subplot.bottom'] = 0.145
        self.parameters['subplot.right'] = 0.9
        
        super().__init__(**self.parameters )
        

#---------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------- 
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------- 
class TexPaper(StandardImproved):

    def __init__(self ,**kwargs):
        
        self.parameters = kwargs
      
        kwarg = {'text':'True', 'box':'gray', 'axeslabel':'#555555', 'fontsize':'16', 'legendfontsize':'12' , 'labelsize':[16,14,14]}
        self.parameters.update(kwarg)
    
        #self.parameters['grid.linestyle'] = '--' 
        #self.parameters['grid.alpha'] = '1'
        self.parameters['box'] = 'gray' 
        self.parameters['axeslabel'] = '#2E3436'
        #self.parameters['xtickcolor'] = 'gray' 
        #self.parameters['ytickcolor'] = 'gray' 
        #self.parameters['linestyle'] = 'paper'
        #if 'scheme' not in self.parameters.keys():
        #    self.parameters['scheme'] = 'grayscale'
        self.parameters['legend.loc'] = 'upper right'
        self.parameters['label'] = ['$t$ [s]', '$y(t)= f(x)$']
        self.parameters['linestyle'] = 'ls3a'
        self.parameters['scheme'] = '2gray'
        self.parameters['cycle']  = self.cycle('1')
        self.parameters['subplot.bottom'] = 0.145
        self.parameters['axes.linewidth'] = 0.9
        super().__init__(**self.parameters )


#----------------------------------------------------------------------------------------------------------

class Helvet(StandardImproved):
     
    
     def __init__(self ,**kwargs):
        
        self.parameters = kwargs
     
         

        kwarg = {'box':'gray', 'axeslabel':'#2E3436', 'fontsize':'14', 'legendfontsize':'12' , 'labelsize':[16,14,14], **kwargs }
        self.parameters.update(kwarg)
        
        self.parameters['fontstyle'] = 'normal'
        self.parameters['font'] = 'Helvetica'
        #self.parameters['box'] = 'gray' 
        #self.parameters['axeslabel'] = 'gray'
        #self.parameters['xtickcolor'] = 'gray' 
        #self.parameters['ytickcolor'] = 'gray' 
        self.parameters['fontsize'] = 12.0
        #self.parameters['legendfontsize'] =12.0
        self.parameters['xtick.labelsize'] = 12
        self.parameters['ytick.labelsize'] = 12
        if 'scheme' not in self.parameters.keys():
            self.parameters['scheme'] = 'tthmod'
        #self.parameters['gridcolor'] = 'gray'
        #self.parameters['legend.loc'] = 'upper right'
        self.parameters['lines.linewidth'] = 3
        #self.parameters['grid.linestyle'] = '--'
        #self.parameters['grid.linewidth'] = 0.7
        #self.parameters['grid.dashes'] = (4,4)
        #self.parameters['title.colors'] = '#2E3436'
        #self.parameters['axes.linewidth'] = 0.7
        #super().get({'title.colors':'#2E3436'})
        
        super().__init__(**self.parameters )
        

#---------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------

class HelvetBeamer(StandardImproved):
     
     def __init__(self , **kwargs):
        
        self.parameters = kwargs

        self.parameters.update(kwargs)
        
        self.parameters['fontstyle'] = 'normal'
        self.parameters['font'] = 'Helvetica'
        #self.parameters['box'] = 'gray'    #'gray' 
        self.parameters['axeslabel'] = '#2E3436'
        #self.parameters['xtickcolor'] = 'gray' 
        #self.parameters['ytickcolor'] = 'gray' 
        self.parameters['xtick.labelsize'] = 12
        self.parameters['ytick.labelsize'] = 12
        self.parameters['fontsize'] = 18.0
        self.parameters['legendfontsize'] =14.0
        #self.parameters['xtick.labelsize'] = 26
        #self.parameters['ytick.labelsize'] = 26
        if 'scheme' not in self.parameters.keys():
            self.parameters['scheme'] = 'mycolor'
        #self.parameters['gridcolor'] = 'gray'
        #self.parameters['legend.loc'] = 'best'
        self.parameters['lines.linewidth'] = 4
        #self.parameters['grid.linestyle'] = '--'
        self.parameters['grid.linewidth'] = 0.7
        #self.parameters['grid.dashes'] = (4,4)
        self.parameters['title.colors'] = 'k'
        self.parameters['figuresize'] = [5,4]
        #super().get({'title.colors':'#2E3436'})
        self.parameters['subplot.bottom'] = 0.145
        self.parameters['axes.linewidth'] = 1
        
        super().__init__(**self.parameters )
        

#---------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------

class Times(StandardImproved):
   
     def __init__(self, **kwargs):
        
        self.parameters = kwargs
     

        self.parameters['fontstyle'] = 'normal'
        self.parameters['font.weight'] = 'normal'
        self.parameters['font'] = 'Times New Roman '
        
        if 'fontsize' not in self.parameters.keys():
            self.parameters['fontsize'] = 14.0
            
        if 'legendfontsize' not in self.parameters.keys():
            self.parameters['legendfontsize'] =12.0
        
        if 'xtick.labelsize' not in self.parameters.keys():
            self.parameters['xtick.labelsize'] = 16
        
        if 'ytick.labelsize' not in self.parameters.keys():
            self.parameters['ytick.labelsize'] = 16
         
        if 'box' not in self.parameters.keys():
            self.parameters['box'] = '#AAAAAA' 
        
        if 'lines.linewidth' not in self.parameters.keys():
            self.parameters['lines.linewidth'] = 3
        
        if 'title.colors' not in self.parameters.keys():
            self.parameters['title.colors'] = '#2E3436'
        
        if 'scheme' not in self.parameters.keys():
            self.parameters['scheme'] = 'nb'
        
        if 'legend.loc' not in self.parameters.keys():
            self.parameters['legend.loc'] = 'upper right'
        
        if 'legend.linewidth' not in self.parameters.keys():
            self.parameters['lines.linewidth'] = 2

        if 'legend.framelinewidth' not in self.parameters.keys():
            self.parameters['legend.framelinewidth'] = 0.2
            
        if 'legend.textcolor' not in self.parameters.keys():
            self.parameters['legend.textcolor'] = '#2E3436'
       
        super().__init__(**self.parameters )

#---------------------------------------------------------------------------------------------------------

class TimesItalic(StandardImproved):
   
     def __init__(self, **kwargs):
        
        self.parameters = kwargs

        
        fontname = os.path.join(matplotlib.get_data_path(),
                            'fonts', 'ttf', 'TimesNewRomanItalic.ttf')

        self.parameters.update(kwargs)
        
        self.parameters['fontstyle'] = 'italic'
        self.parameters['font'] = 'Times New Roman'
        self.parameters['font.sans-serif'] = 'Times New Roman' 
        self.parameters['box'] = '#AAAAAA' 
        
        super().__init__(**self.parameters )

#---------------------------------------------------------------------------------------------------------



class SansItalic(StandardImproved):
   
    def __init__(self, **kwargs):
        
        self.parameters = kwargs

        self.parameters.update(kwargs)
        self.parameters['fontstyle'] = 'italic'
        self.parameters['font'] = 'serif'
        self.parameters['font.sans-serif'] = 'sans-serif' 
        self.parameters['box'] = 'gray' 
        self.parameters['axeslabel'] = 'gray'
        self.parameters['xtickcolor'] = 'gray' 
        self.parameters['ytickcolor'] = 'gray' 
        self.parameters['fontsize'] = 12.0
        self.parameters['legendfontsize'] =12.0
        self.parameters['xtick.labelsize'] = 16
        self.parameters['ytick.labelsize'] = 16
        self.parameters['scheme'] = 'nb'
        self.parameters['gridcolor'] = 'gray'   # '#dddddd' # 'gray'  #  #  #
        self.parameters['legend.loc'] = 'upper left'
        self.parameters['lines.linewidth'] = 3
        self.parameters['grid.dashes'] = (5,5)       
        self.parameters['title.colors'] = '#2E3436'
        #self.parameters['legend.bbox'] = (1,0.5)
        self.parameters['legend.loc'] = 'best' 
        if 'scheme' not in self.parameters.keys():
            self.parameters['scheme'] = 'nb2'   #'nb2' #'grayimproved' #'nb2'
        self.parameters['grid.linewidth'] = 0.7
        #self.parameters['legend.axespad'] = 0. 
        self.parameters['axes.linewidth'] = 0.7

        super().__init__(**self.parameters )
    
       # #82A6C7 #DAADEB #99D366

#---------------------------------------------------------------------------------------------------------







class PalatinoItalic(StandardImproved):
   
     def __init__(self ,**kwargs):
        
        self.parameters = kwargs
        self.parameters.update(kwargs)
        
        self.parameters['fontstyle'] = 'italic'
        self.parameters['font'] = 'Palatino-Italic'
        
        self.parameters['axeslabel'] = '#2E3436'
        self.parameters['xtickcolor'] = 'gray' 
        self.parameters['ytickcolor'] = 'gray' 
        
        if 'fontsize' not in self.parameters.keys():
            self.parameters['fontsize'] = 18.0
        if 'legendfontsize' not in self.parameters.keys():
            self.parameters['legendfontsize'] =14.0
        if 'xtick.labelsize' not in self.parameters.keys():
            self.parameters['xtick.labelsize'] = 14
        if 'ytick.labelsize' not in self.parameters.keys():
            self.parameters['ytick.labelsize'] = 14
        if 'scheme' not in self.parameters.keys():
            self.parameters['scheme'] = 'nb'
            self.parameters['gridcolor'] = 'gray'
            self.parameters['legend.loc'] = 'best'
        if 'lines.linewidth' not in self.parameters.keys():
            self.parameters['lines.linewidth'] = 3
        if 'grid.dashes' not in self.parameters.keys():       
            self.parameters['grid.dashes'] = (8,8)       
        
        
        self.parameters['subplot.bottom'] = 0.145
        self.parameters['grid.linewidth'] = 0.7
        self.parameters['title.colors'] = '#2E3436'
        self.parameters['axes.linewidth'] = 0.7
        self.parameters['legendEdgeColor'] = 'gray'
        
        super().__init__(**self.parameters )


#---------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------------


class PalatinoItalicLight(BasePlot):
   
     def __init__(self ,**kwargs):
        
        self.parameters = kwargs
        
        self.parameters.update(kwargs)
        self.parameters['fontstyle'] = 'italic'
        self.parameters['font'] = 'Palatino-Italic'
        #self.parameters['font.sans-serif'] = 'Times New Roman' 
        self.parameters['box'] = 'gray' 
        self.parameters['axeslabel'] = '#2E3436'
        self.parameters['xtickcolor'] = 'gray' 
        self.parameters['ytickcolor'] = 'gray' 
        if 'fontsize' not in self.parameters.keys():
            self.parameters['fontsize'] = 18.0
        if 'legendfontsize' not in self.parameters.keys():
            self.parameters['legendfontsize'] =14.0
        if 'xtick.labelsize' not in self.parameters.keys():
            self.parameters['xtick.labelsize'] = 14
        if 'ytick.labelsize' not in self.parameters.keys():
            self.parameters['ytick.labelsize'] = 14
        if 'scheme' not in self.parameters.keys():
            self.parameters['scheme'] = 'nb'
        self.parameters['gridcolor'] = 'gray'
        self.parameters['legend.loc'] = 'best'
        if 'lines.linewidth' not in self.parameters.keys():
            self.parameters['lines.linewidth'] = 2
        self.parameters['grid.dashes'] = (6,6)       
        self.parameters['subplot.bottom'] = 0.145
        self.parameters['grid.linewidth'] = 0.5
        self.parameters['title.colors'] = '#2E3436'
        self.parameters['axes.linewidth'] = 0.5
        self.parameters['legendEdgeColor'] = 'gray'
        
        super().__init__(**self.parameters )


#---------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------


class PalatinoItalicDark(BasePlot):
   
     def __init__(self ,**kwargs):
        
        self.parameters = kwargs
        
        self.parameters.update(kwargs)
        self.parameters['fontstyle'] = 'italic'
        self.parameters['font'] = 'Palatino-Italic'
        #self.parameters['font.sans-serif'] = 'Times New Roman' 
        self.parameters['box'] = 'gray' 
        self.parameters['axeslabel'] = '#2E3436'
        self.parameters['xtickcolor'] = 'gray' 
        self.parameters['ytickcolor'] = 'gray' 
        if 'fontsize' not in self.parameters.keys():
            self.parameters['fontsize'] = 18.0
        if 'legendfontsize' not in self.parameters.keys():
            self.parameters['legendfontsize'] =14.0
        if 'xtick.labelsize' not in self.parameters.keys():
            self.parameters['xtick.labelsize'] = 14
        if 'ytick.labelsize' not in self.parameters.keys():
            self.parameters['ytick.labelsize'] = 14
        if 'scheme' not in self.parameters.keys():
            self.parameters['scheme'] = 'nb'
        self.parameters['gridcolor'] = 'gray'
        self.parameters['legend.loc'] = 'best'
        if 'lines.linewidth' not in self.parameters.keys():
            self.parameters['lines.linewidth'] = 2
        self.parameters['grid.dashes'] = (4,4)       
        self.parameters['subplot.bottom'] = 0.145
        self.parameters['grid.linewidth'] = 1
        self.parameters['title.colors'] = '#2E3436'
        self.parameters['axes.linewidth'] = 1
        self.parameters['legendEdgeColor'] = 'gray'
        self.parameters['tick.MajW'] = 1
        self.parameters['tick.minW'] = 1

        super().__init__(**self.parameters )


#---------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------


class Palatino(BasePlot):
     
     def __init__(self , **kwargs):
        
        self.parameters = kwargs
        
        self.parameters['fontstyle'] = 'italic'
        self.parameters['font'] = 'Palatino'
        self.parameters['box'] = 'gray' 
        self.parameters['axeslabel'] = 'gray'
        self.parameters['xtickcolor'] = 'gray' 
        self.parameters['ytickcolor'] = 'gray' 
        self.parameters['fontsize'] = 14.0
        self.parameters['legendfontsize'] =12.0
        self.parameters['xtick.labelsize'] = 12
        self.parameters['ytick.labelsize'] = 12
        if 'scheme' not in self.parameters.keys():
            self.parameters['scheme'] = 'nb'
        self.parameters['gridcolor'] = 'gray'
        self.parameters['legend.loc'] = 'best'
        self.parameters['lines.linewidth'] = 4
        self.parameters['grid.dashes'] = (4,4)       
        self.parameters['grid.linewidth'] = 0.7
        self.parameters['axes.linewidth'] = 0.7
        
        super().__init__(**self.parameters )


#---------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------

class PalatinoImproved(BasePlot):
     
     def __init__(self , **kwargs):
        
        self.parameters = kwargs
        
        self.parameters['fontstyle'] = 'italic'
        self.parameters['font'] = 'Palatino-italic'
        self.parameters['box'] = 'gray' 
        self.parameters['axeslabel'] = '#2E3436'
        self.parameters['xtickcolor'] = 'gray' 
        self.parameters['ytickcolor'] = 'gray' 
        self.parameters['fontsize'] = 14.0
        self.parameters['legendfontsize'] =12.0
        self.parameters['xtick.labelsize'] = 12
        self.parameters['ytick.labelsize'] = 12
        if 'scheme' not in self.parameters.keys():
            self.parameters['scheme'] = 'nb'
        self.parameters['gridcolor'] = 'gray'
        self.parameters['legend.loc'] = 'best'
        self.parameters['lines.linewidth'] = 3
        self.parameters['grid.dashes'] = (3,3)       
        self.parameters['title.colors'] = '#2E3436'
        self.parameters['grid.linewidth'] = 1
        self.parameters['axes.linewidth'] = 0.7
        
        super().__init__(**self.parameters )


#---------------------------------------------------------------------------------------------------------


class DustimoFont(BasePlot):
     
     def __init__(self , **kwargs):
        
        self.parameters = kwargs
        
        self.parameters['fontstyle'] = 'italic'
        #self.parameters['font'] = 'RomanSerif'
        #self.parameters['font'] = 'Thryomanes'
        self.parameters['font'] = 'Dustismo Roman'
        #self.parameters['font'] = 'Brougham'
        #self.parameters['font.sans-serif'] = 'Times New Roman' 
        self.parameters['box'] = 'gray' 
        self.parameters['axeslabel'] = 'gray'
        self.parameters['xtickcolor'] = 'gray' 
        self.parameters['ytickcolor'] = 'gray' 
        self.parameters['fontsize'] = 12.0
        self.parameters['legendfontsize'] =12.0
        self.parameters['xtick.labelsize'] = 10
        self.parameters['ytick.labelsize'] = 10
        if 'scheme' not in self.parameters.keys():
            self.parameters['scheme'] = 'nb'
        self.parameters['gridcolor'] = 'gray'
        self.parameters['legend.loc'] = 'upper right'
        self.parameters['lines.linewidth'] = 3
        self.parameters['grid.dashes'] = (3,3)       
        self.parameters['axes.linewidth'] = 1
        self.parameters['grid.linewidth'] = 1
        
        super().__init__(**self.parameters )



#---------------------------------------------------------------------------------------------------------

class Fonts(BasePlot):
     
     def __init__(self , **kwargs):
        
        self.parameters = kwargs
        
        #self.parameters['fontstyle'] = 'italic'
        self.parameters['fontstyle'] = 'italic'
        self.parameters['font'] = 'Times Sans Serif'
        self.parameters['font.sans-serif'] = 'Times Sans Serif' 
        #self.parameters['font'] = 'RomanSerif'
        #self.parameters['font'] = 'Thryomanes'
        #self.parameters['font'] = 'Dustismo Roman'
        #self.parameters['font'] = 'Brougham'
        #self.parameters['font'] = 'Kepler Std'
        #self.parameters['font.sans-serif'] = 'Times New Roman' 
        #self.parameters['font'] = 'Palatino-italic'
        self.parameters['box'] = 'gray' 
        self.parameters['axeslabel'] = '#2E3436'
        self.parameters['xtickcolor'] = 'gray' 
        self.parameters['ytickcolor'] = 'gray' 
        self.parameters['fontsize'] = 14.0
        self.parameters['legendfontsize'] =12.0
        self.parameters['xtick.labelsize'] = 16
        self.parameters['ytick.labelsize'] = 16
        if 'scheme' not in self.parameters.keys():
            self.parameters['scheme'] = 'nb'
        self.parameters['gridcolor'] = 'gray'
        self.parameters['legend.loc'] = 'upper right'
        self.parameters['lines.linewidth'] = 3
        self.parameters['grid.dashes'] = (6,6)       
        self.parameters['axes.linewidth'] = 0.7
        self.parameters['grid.linewidth'] = 0.7
        self.parameters['title.colors'] = '#2E3436'
        
        super().__init__(**self.parameters )

##################################################################
##################################################################
##################################################################

class FontsDark(BasePlot):
     
     def __init__(self , **kwargs):
        
        self.parameters = kwargs
        
        #self.parameters['fontstyle'] = 'italic'
        self.parameters['fontstyle'] = 'italic'
        self.parameters['font'] = 'Times Sans Serif'
        self.parameters['font.sans-serif'] = 'Times Sans Serif' 
        #self.parameters['font'] = 'RomanSerif'
        #self.parameters['font'] = 'Thryomanes'
        #self.parameters['font'] = 'Dustismo Roman'
        #self.parameters['font'] = 'Brougham'
        #self.parameters['font'] = 'Kepler Std'
        #self.parameters['font.sans-serif'] = 'Times New Roman' 
        #self.parameters['font'] = 'Palatino-italic'
        self.parameters['box'] = 'gray' 
        self.parameters['axeslabel'] = '#2E3436'
        self.parameters['xtickcolor'] = 'gray' 
        self.parameters['ytickcolor'] = 'gray' 
        self.parameters['fontsize'] = 14.0
        self.parameters['legendfontsize'] =12.0
        self.parameters['xtick.labelsize'] = 16
        self.parameters['ytick.labelsize'] = 16
        if 'scheme' not in self.parameters.keys():
            self.parameters['scheme'] = 'nb'
        self.parameters['gridcolor'] = 'gray'
        self.parameters['legend.loc'] = 'upper right'
        self.parameters['lines.linewidth'] = 3
        self.parameters['grid.dashes'] = (3,3)       
        self.parameters['axes.linewidth'] = 1
        self.parameters['grid.linewidth'] = 1
        self.parameters['title.colors'] = '#2E3436'
        
        super().__init__(**self.parameters )




             


