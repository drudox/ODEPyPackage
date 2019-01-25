# -*- coding : utf-8 -*- 

from abc import ABCMeta, abstractmethod
import matplotlib.pyplot as plt
import numpy as np
from abc import ABCMeta, abstractmethod
from cycler import cycler
import matplotlib.font_manager as font_manager
from matplotlib.ticker import MultipleLocator, FormatStrFormatter 
from matplotlib.axes import Axes
import matplotlib
from matplotlib import cm
from collections import OrderedDict
from matplotlib.ticker import AutoMinorLocator


class BasePlot(metaclass=ABCMeta):
   
   @abstractmethod  
   def __init__(self,**kwargs):
      self.parameters = kwargs   
      self.nrows = None
      self.ncols = None
      self.figsize = None
      
      font_dirs = ['/home/marco/.fonts', ]
      font_files = font_manager.findSystemFonts(fontpaths=font_dirs)
      font_list = font_manager.createFontList(font_files)
      font_manager.fontManager.ttflist.extend(font_list)
      
   
   #----------------------------------------------------------------------------------------------
   def findLimits(self):
         '''
             found the vectors limits in order to set the figure limits x,y
         '''
         self._minX = 10E20
         self._minY = 10E20
         self._maxX = 0.000
         self._maxY = 0.000
         for i in range(0,len(self.variable)-1,3):
            if self._minX >= min(self.variable[i]):  
               self._minX = min(self.variable[i])
            if self._maxX <= max(self.variable[i]):
               self._maxX = max(self.variable[i])
            if self._minY >= min(self.variable[i+1]):  
               self._minY = min(self.variable[i+1])
            if self._maxY <= max(self.variable[i+1]):
               self._maxY = max(self.variable[i+1])
         return [round(self._minX,2), round(self._maxX,2) , self._minY, self._maxY]
         #return [round(self._minX,2), round(self._maxX,2) , round(self._minY,2) , round(self._maxY,2)]

      
   def schemes(self, style:str = 'nb'):
         
       if style == 'nb':
            return ['#8DA0CB', '#E58AC3', '#A6D853', '#FFD930', '#B3B3B3', '#5FC3A4', '#FC8D62', '#E5C494']
       elif style == 'nb2': ##7aa0c4","#ca82e1" ,"#8bcd50
            return  ['#7aa0c4', '#ca82e1', '#8bcd50', '#FFD930', '#B3B3B3', '#5FC3A4', '#FC8D62', '#66C2A5'] 
       elif style == 'nbc': ##7aa0c4","#ca82e1" ,"#8bcd50
            return  ['#7aa0c4', '#ca82e1', '#8bcd50','#B3B3B3', '#FFD930', '#5FC3A4', '#FC8D62', '#66C2A5'] 
       elif style == 'nbd':
            return ['#8DA0CB', '#E58AC3', '#A6D853', '#B3B3B3','#FFD930', '#5FC3A4', '#FC8D62', '#E5C494']
       elif style == 'mycolor': ##7aa0c4","#ca82e1" ,"#8bcd50
            return  ['#7aa0c4', '#ca82e1', '#8bcd50',"#FFD92F", '#B3B3B3',"#64b9a1", "#745ea6", "#db7e76" ]
       elif style == 'migcolor': 
            return  ['#66C2A5', '#FC8D62', '#8DA0CB', '#E78AC3', '#A6D854', '#FFD92F', '#E5C494', '#B3B3B3']
       elif style == 'vega':
            return  ['#1F77B4', '#FF7F0E', '#2CA02C', '#D62728', '#9467BD', '#8C564B', '#E377C2', '#7F7F7F', '#BCBD22', '#17BECF']
       elif style ==  'gg':
            return ['#E24A33', '#348ABD', '#988ED5', '#777777', '#FBC15E', '#8EBA42', '#FFB5B8']     
       elif style == 'brewers':  
            return  ['#66C2A5', '#FC8D62', '#8DA0CB', '#E78AC3', '#A6D854', '#FFD92F', '#E5C494', '#B3B3B3']
       elif style == 'nb2':
            return  ['#8DA0CB', '#E78AC3', '#A6D854', '#FFD92F', '#E5C494', '#B3B3B3', '#66C2A5', '#FC8D62']
        # color['nb3'] =  ['#8DA0CB', '#E78AC3', '#A6D854', '#FFD92F', '#E5C494', '#B3B3B3', '#66C2A5', '#FC8D62']
       elif style == 'mystyle':
            return  ['#8DA0CB', '#E58AC3', '#A6D853', '#FFD930', '#B2B2B2', '#5FC3A4', '#FC8D62', '#66C2A5']
       elif style == 'tthmod':
            return  ['#30a2da', '#fc4f30', '#e5ae38', '#6d904f', '#8b8b8b']
       elif style == 'grayscale':
            return  ['#252525', '#525252', '#737373', '#969696', '#BDBDBD', '#D9D9D9', '#F0F0F0', '#F0F0FF' ]
       elif style == 'grayscale2':
            return  ['#525252', '#969696', '#BDBDBD' ,'#D9D9D9', '#F0F0F0', '#F0F0FF' ]
       elif style == 'grayscale3':
            return  ['#4B4B4B', '#686868', '#868686' ,'#B2B2B2', '#CFCFCF', '#F0F0FF' ]
       elif style =='accent':
            return  ['#7fc97f','#beaed4','#fdc086','#ffff99','#386cb0','#f0027f','#bf5b17','#666666']
       elif style == 'palette2':
            return  ['#1B9E77','#D95F02','#7570B3','#E7298A','#66A61E','#E6AB02','#A6761D','#666666']
       elif style == 'set1':
            return  ['#E41A1C','#377EB8','#4DAF4A','#984EA3','#FF7F00','#FFFF33', '#F781BF','#A65628']
       elif style == 'set2':
            return  ['#66C2A5','#FC8D62','#8DA0CB','#E78AC3','#A6D854','#FFD92F','#E5C494','#B3B3B3']
       elif style =='set3':
            return  ['#8dd3c7','#ffd92f','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5']
       elif style =='pastel1':
            return ['#b3e2cd','#fdcdac','#cbd5e8','#f4cae4','#e6f5c9','#fff2ae','#f1e2cc','#cccccc']
       elif style == 'pastel2':
            return ['#fbb4ae','#b3cde3','#ccebc5','#decbe4','#fed9a6','#ffffcc','#e5d8bd','#fddaec']
       elif style == 'paraview':
            return ['#3B4CC0','#7C9FF9','#C0D4F5','#F2CBB7','#EE8568','#B40426']
       elif style =='paraview2':
            return ['#3B4CC0','#6889EE','#9ABAFF','#C9D8EF','#EDD1C2','#F7A889','#E26A53','#B40426']
       elif style =='radiance1':
            return  ['#EDED49','#E6B500','#E66300','#E61000','#A40000','#520000','#000000']
       elif style =='xray':
            return  ['#000000','#333333','#666666','#999999','#CCCCCC']
       elif style == 'bluesea':
            return  ['#22204D','#35426B','#3C778F','#4D9BA8','#5DBCBF','#70D4CD','#9DEADA']
       elif style == 'bluesky':
            return  ['#5C5C5C','#59647F','#4D67A0','#4872C4','#5085E3','#62A1F7','#86C1FF','#BFE0FF']
       elif style =='bluenice':
            return  ['#3B4CC0','#4F69D9','#6585EC','#7C9FF9','#93B5FF','#AAC7FD','#C0D4F5','#D4DBE6']
       elif style =='paired':
            return  ['#A6CEE3','#1F78B4','#B2DF8A','#33A02C','#FB9A99','#E31A1C','#FDBF6F','#FF7F00','#CAB2D6','#6A3D9A']
       elif style == 'new':
            return  ['#30A2DA','#FC4F30','#6D904F','#E5AE38','#8B8B8B','#FFB5B8','#988ED5','#30A2DA']
       elif style == 'new2':
            return  ['#0072B2','#A60628','#467821','#F0E442','#7A68A6','#56B4E9','#CC79A7','#009E73']
       elif style == 'new3':
            return  ['#0072B2','#A60628','#467821','#F0E442','#7A68A6','#56B4E9','#CC79A7','#009E73','#348ABD','#D55E00']
       elif style == 'new4':
            return  ['#92C6FF','#FF9F9A','#97F0AA','#FFD92F','#B0E0E6','#D0BBFF','#B8BEAF','#FFB6FC']
       elif style == 'new5':
            return  ['#7A8FB1','#EC9CAE','#B2E295','#FFDFA8','#405C89','#BE4F69','#72B54B','#CFA356']
       elif style == 'new6':
            return  ['#B8CCED','#FBBDCB','#D2F8BB','#FFDA9A']
       elif style =='new6a':
            return  ['#B8CCED','#FBBDCB']
       elif style == '3gray':
            return ['#000000' , '#000000', '#000000' , '#666666','#666666','#666666', '#999999', '#999999', '#999999']      
       elif style == '2gray':
            return ['#000000' ,'#888A85', '#BABDB6', '#000000','#888A85', '#BABDB6']     
       elif style == 'grayimproved':
            return ['#000000','#808080','#778899','#A9A9A9','#B0C4DE','#C0C0C0']
       elif style == 'monocrome':
            return ['#000000','#666666','#999999','#B3B3B3']

#------------------------------------------------------------------------------------------
   
   def linestyles(self, style : str = 'ls1'):           
         
       if style =='ls1':
           return  ['-', '--', ':', '-.']
       elif style == 'ls2':             
           return  ['-', '--', ':', '-.','-','--',':','-.'] 
       elif style == 'ls3':
           return  ['-','-','-','-.','-.','-.',':',':',':']       
       elif style == 'ls3a':
           return  ['-','-','-.','-.',':',':']       
       elif style == 'ls10':
           return ['-', '--', ':', '-.','-','--',':','-.','-','--'] 
       elif style =='lsa8':
           return ['-','--','-','--','-','--','-','--']
       elif style == 'lp8':
           return ['-',':','-',':','-',':','-',':']
       elif style == 'llt10':
           return ['-','--','-','--','-','--','-','--','-','--']
       elif style == 'lp10':
           return ['-',':','-',':','-',':','-',':','-',':']
       elif style == 'ldash':
           return [(0, ()),(0, (1, 10)),(0, (1, 5)),(0, (1, 1))]
       elif style == 'ls8':
           return [(0, ()),(0, (1, 1)),(0, (5, 10)),(0, (5, 5)),(0, (5, 1)),(0, (3, 10, 1, 10)),(0, (3, 5, 1, 5)),(0, (3, 5, 1, 5, 1, 5)) ]
       elif style == 'ln':
           return ["--", (0,(5,2,5,5,1,4))]
       elif style == 'paper':
           return [(0, ()), (0, (3, 1, 1, 1, 1, 1)), (0, (5, 5)),(0, (5, 1)),(0, (3, 5, 1, 5)),(0, (3, 5, 1, 5, 1, 5)),(0, (1, 1)),(0, (5, 10)) ]
         
         
        # linestyle = OrderedDict(
        #    [('solid',               (0, ())),
        #    ('loosely dotted',      (0, (1, 10))),
        #    ('dotted',              (0, (1, 5))),
        #    ('densely dotted',      (0, (1, 1))),
        #    ('loosely dashed',      (0, (5, 10))),
        #    ('dashed',              (0, (5, 5))),
        #    ('densely dashed',      (0, (5, 1))),
        #    ('loosely dashdotted',  (0, (3, 10, 1, 10))),
        #    ('dashdotted',          (0, (3, 5, 1, 5))),
        #    ('densely dashdotted',  (0, (3, 1, 1, 1))),
        #    ('loosely dashdotdotted', (0, (3, 10, 1, 10, 1, 10))),
        #    ('dashdotdotted',         (0, (3, 5, 1, 5, 1, 5))),
        #    ('densely dashdotdotted', (0, (3, 1, 1, 1, 1, 1)))])
#----------------------------------------------------------------------------------------------------
   def cycle(self,n : str):
            if n == '0':
                  return plt.cycler("color", self.schemes(self.parameters['scheme']) )      #colors)         
            elif n=='1':
                  return plt.cycler("color", self.schemes(self.parameters['scheme']) ) + plt.cycler("linestyle", self.linestyles(self.parameters['linestyle']))
            elif n=='2':
                  return plt.cycler("color", self.schemes(self.parameters['scheme']) ) * plt.cycler("linestyle", self.linestyles(self.parameters['linestyle']))
            elif n=='3':
                  return plt.cycler("linestyle", self.linestyles(self.parameters['linestyle'])) * plt.cycler("color", self.schemes(self.parameters['scheme']))
            else :   
                 raise ValueError('Index out of bound')





#------------------------------------------------------------------------------------------------
   def setparams(self):
             
             #self.parameters.update(kwargs)
             
             if 'box' not in self.parameters.keys() :
                self.parameters['box'] = '#AAAAAA' 
             if 'axeslabel' not in self.parameters.keys() :
                self.parameters['axeslabel'] = 'gray'
             if 'xtickcolor' not in self.parameters.keys() :
                self.parameters['xtickcolor'] = 'gray' 
             if 'ytickcolor' not in self.parameters.keys() :
                self.parameters['ytickcolor'] = 'gray' 
             if 'gridcolor' not in self.parameters.keys():
                self.parameters['gridcolor'] = '#dddddd'
             if 'font' not in self.parameters.keys():
                self.parameters['font'] = 'serif'
             if 'fontstyle' not in self.parameters.keys():
                self.parameters['fontstyle'] = 'italic'
             if 'text' not in self.parameters.keys():
                self.parameters['text'] = False           
             if 'fontsize' not in self.parameters.keys():   
                self.parameters['fontsize'] = 10.0
             if 'legendfontsize' not in self.parameters.keys():   
                self.parameters['legendfontsize'] =10.0
             if 'legendEdgeColor' not in self.parameters.keys():
                self.parameters['legendEdgeColor'] = '#dddddd'
             if 'labelsize' not in self.parameters.keys(): 
                self.parameters['labelsize'] = [12,12,12]
             if 'linestyle' not in self.parameters.keys():
                self.parameters['linestyle'] = self.linestyles('paper')
             if 'scheme' not in self.parameters.keys():
                self.parameters['scheme'] = 'nb'
             if 'cycle' not in self.parameters.keys():
                self.parameters['cycle']  = self.cycle('0')
             if 'legend.loc' not in self.parameters.keys():
                self.parameters['legend.loc'] = 'best'
             if 'lines.linewidth' not in self.parameters.keys():
                self.parameters['lines.linewidth'] = '2' 
             if 'grid.alpha' not in self.parameters.keys():
                self.parameters['grid.alpha'] = '1'
             if 'text.color' not in self.parameters.keys():
                self.parameters['text.color'] = 'gray'
             if 'font.sans-serif' not in self.parameters.keys():   
                self.parameters['font.sans-serif'] = 'serif'
             if  'figure.figsize' not in self.parameters.keys():
                self.parameters['figure.figsize'] = [9.5,4.5]
             if 'subplot.bottom' not in self.parameters.keys():
                self.parameters['subplot.bottom'] = 0.125
             if 'text' not in self.parameters.keys():
                self.parameters['text'] = False
             if 'package' not in self.parameters.keys():   
                self.parameters['package'] = []
             if 'grid.linewidth' not in self.parameters.keys():
                self.parameters['grid.linewidth'] = 1
             if 'axes.linewidth' not in self.parameters.keys(): 
                self.parameters['axes.linewidth'] = 0.7
             if 'font.weight' not in self.parameters.keys():
                self.parameters['font.weight'] = 'normal'
             if 'legend.borderpad' not in self.parameters.keys():   
                self.parameters['legend.borderpad'] = 0.5
             if 'legend.borderaxespad' not in self.parameters.keys():
                self.parameters['legend.borderaxespad'] =0.5
             if 'legend.framelinewidth' not in self.parameters.keys():  
                self.parameters['legend.framelinewidth'] = 1.0
             if 'boxplot.boxprops.linestyle' not in self.parameters.keys():
                self.parameters['boxplot.boxprops.linestyle'] = '--'
             if 'tick.MajW' not in self.parameters.keys(): 
                self.parameters['tick.MajW'] = 0.7
             if 'tick.minW' not in self.parameters.keys(): 
                self.parameters['tick.minW'] = 0.7
                       
    




             myparams = {
               'patch.linewidth' : '0.5',
               'patch.facecolor' : '#348ABD',  # blue
               'patch.edgecolor' : 'EEEEEE',
               'patch.antialiased' : True,
               'font.size': self.parameters['fontsize'],
               'xtick.major.size' : 0.1,
               'xtick.minor.size' : 0.025,
               'axes.edgecolor': self.parameters['box'] ,   
               'axes.linewidth': self.parameters['axes.linewidth'],   # BOX width
               'axes.xmargin': 0,    
               'axes.ymargin': 0,     
               'axes.labelcolor': self.parameters['axeslabel'],     
               'axes.axisbelow': True,   
               'xtick.color': self.parameters['xtickcolor'],      
               'ytick.color': self.parameters['ytickcolor'],     
               'axes.prop_cycle': self.parameters['cycle'],     
               'grid.linestyle': '--', 
               'grid.alpha': self.parameters['grid.alpha'],
               'grid.linewidth' : self.parameters['grid.linewidth'],
               'grid.color' :  self.parameters['gridcolor'],
               'font.family': self.parameters['font'] ,
               #'labelsize' : self.parameters['labelsize'],
               'font.size'  : self.parameters['fontsize'],
               'legend.edgecolor' : self.parameters['legendEdgeColor'],
               'legend.fancybox'  : False,
               'legend.borderpad' : self.parameters['legend.borderpad'],
               'legend.fontsize'  : self.parameters['legendfontsize'], 
               'legend.loc'       : self.parameters['legend.loc'] ,
               'legend.framealpha': 1,
               'legend.borderaxespad' : self.parameters['legend.borderaxespad'],
               'lines.linewidth' : self.parameters['lines.linewidth'],
               'font.style' : self.parameters['fontstyle'],
               'figure.autolayout': False,
               'font.sans-serif' : self.parameters['font.sans-serif'], 
               'text.usetex': self.parameters['text'] ,
               'patch.linewidth' : self.parameters['legend.framelinewidth'],
               'boxplot.boxprops.linestyle' : self.parameters['boxplot.boxprops.linestyle'],
               #'figure.dpi': 72.0,
               #'figure.edgecolor': (1, 1, 1, 0),
               #'figure.facecolor': (1, 1, 1, 0),
               'figure.figsize': self.parameters['figure.figsize'], #[6.0, 4.0],
               #'figure.frameon': True,
               #'figure.max_open_warning': 20,
               'figure.subplot.bottom': self.parameters['subplot.bottom'], #0.125,
               'figure.subplot.hspace': 0.2,
               'figure.subplot.left': 0.125,
               'figure.subplot.right': 0.9,
               'figure.subplot.top': 0.9,
               'figure.subplot.wspace': 0.2,
               'figure.titlesize': 'large',
               'figure.titleweight': 'normal', 
               'text.color' : self.parameters['text.color'],
               #'axes.axisbelow': False,
               #'axes.edgecolor': 'k',
               #'axes.facecolor': 'w',
               #'axes.formatter.limits': [-7, 7],
               #'axes.formatter.use_locale': False,
               #'axes.formatter.use_mathtext': False,
               #'axes.formatter.useoffset': True,
               #'axes.grid': False,
               'axes.grid.axis': 'both',
               #'axes.grid.which': 'major',
               #'axes.hold': True,
               #'axes.labelcolor': 'k',
               #'axes.labelpad': 5.0,
               'axes.labelsize': 'large',
               #'axes.labelweight': 'normal',
               #'axes.linewidth': 1.0,
               #'axes.prop_cycle': cycler('color', ['b', 'g', 'r', 'c', 'm', 'y', 'k']),
               #'axes.spines.bottom': True,
               #'axes.spines.left': True,
               #'axes.spines.right': True,
               #'axes.spines.top': True,
               'axes.titlesize': 'large',
               #'axes.titleweight': 'normal',
               #'axes.unicode_minus': True,
               #'axes.xmargin': 0.0,
               #'axes.ymargin': 0.0,
               #font.family         : sans-serif
               #font.style          : normal
               #font.variant        : normal
               'font.weight'         : self.parameters['font.weight'],
                #font.stretch        : normal
               # 'font.cursive': ['Apple Chancery',
               #            'Textile',
               #            'Zapf Chancery',
               #            'Sand',
               #            'Script MT',
               #            'Felipa',
               #            'cursive'],
               # 'font.family': ['sans-serif'],
               # 'font.fantasy': ['Comic Sans MS',
               #            'Chicago',
               #            'Charcoal',
               #            'ImpactWestern',
               #            'Humor Sans',
               #            'fantasy'],
               # 'font.monospace': ['Bitstream Vera Sans Mono',
               #              'DejaVu Sans Mono',
               #              'Andale Mono',
               #              'Nimbus Mono L',
               #              'Courier New',
               #              'Courier',
               #              'Fixed',
               #              'Terminal',
               #              'monospace'],
               #  'font.sans-serif': ['Bitstream Vera Sans',
               #               'DejaVu Sans',
               #               'Lucida Grande',
               #               'Verdana',
               #               'Geneva',
               #               'Lucid',
               #               'Arial',
               #               'Helvetica',
               #               'Avant Garde',
               #               'sans-serif'],
               #   'font.serif': ['Bitstream Vera Serif',
               #          'DejaVu Serif',
               #          'New Century Schoolbook',
               #          'Century Schoolbook L',
               #          'Utopia',
               #          'ITC Bookman',
               #          'Bookman',
               #          'Nimbus Roman No9 L',
               #          'Times New Roman',
               #          'Times',
               #          'Palatino',
               #          'Charter',
               #          'serif'],
               #    'font.size': 10.0, 
               #     'text.antialiased': True,
               #     'text.color': 'k',
               #     'text.dvipnghack': None,
               #     'text.hinting': 'auto',
               #     'text.hinting_factor': 8,
                    'text.latex.preamble': self.parameters['package'],
               #     'text.latex.preview': False,
               #     'text.latex.unicode': False,
               # 
               # 
               #     'xtick.color': 'k',
               #      'xtick.direction': 'in',
               #      'xtick.labelsize': 'medium',
               #      'xtick.major.pad': 4.0,
               #      'xtick.major.size': 4.0,
               #      'xtick.major.width': 1.5,
               #      'xtick.minor.pad': 4.0,
               #      'xtick.minor.size': 2.0,
               #      'xtick.minor.visible': False,
               #      'xtick.minor.width': 0.5,
               #      'ytick.color': 'k',
               #      'ytick.direction': 'in',
               #      'ytick.labelsize': 'medium',
               #      'ytick.major.pad': 4.0,
               #      'ytick.major.size': 4.0,
               #      'ytick.major.width': 1.5
               #      'ytick.minor.pad': 4.0,
               #      'ytick.minor.size': 2.0,
               #      'ytick.minor.visible': False,
               #      'ytick.minor.width': 0.5
               }
            
            # if self.parameters['text']: 
            #    plt.rc('text', usetex=True )
            #    plt.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
            #    plt.rc('font',family='' ,size=14 )
            #  
            # if 'package' in self.parameters.keys():
            #      for i in self.parameters['package']:
            #         string = r"\usepackage{" + i + "}" 
            #         plt.rcParams['text.latex.preamble']=[string]


             ###########################################
             plt.rcParams.update(myparams)
             ###########################################
    
        
       

#----------------------------------------------------------------------------------------------------------
   def set(self, nrows, ncols, fig , axs, kwargs):
       
       self.config = kwargs
       #print(self.config)

       if 'grid.linestyle' not in self.config.keys():
           self.config['grid.linestyle'] = '--'
       if 'grid.dashes' not in self.config.keys():
            self.config['grid.dashes'] = (5,7)
       if 'title.colors' not in self.config.keys():
            self.config['title.colors'] = '#555555'
       if 'subplot.bottom' not in self.config.keys():
            self.config['subplot.bottom'] = 0.125
       if 'subplot.right' not in self.config.keys():
            self.config['right.bottom'] = 0.93
       if 'legend.bbox' not in self.config.keys():
            self.config['legend.bbox'] = (0,0)
       if 'legend.loc' not in self.config.keys():
            self.config['legend.loc'] = 'best'
       if 'figure.figuresize' not in self.config.keys():
            self.config['figure.figuresize'] = [8,5.5]
       if 'legend.axespad' not in self.config.keys():     
            self.config['legend.axespad'] = 0.
       if 'legend.framelinewidth' not in self.config.keys():    
            self.config['legend.framelinewidth'] = 0.2
       if 'legend.textcolor' not in self.config.keys():
            self.config['legend.textcolor'] = '#555555'
       
       #self.config['legend.axespad'] = 0.

       
       for i in range(nrows):
          for j in range(ncols):
             #self.axs[i,j].set_xlim( self.findLimits()[0] , self.findLimits()[1] )
             #self.axs[i,j].set_ylim( self.findLimits()[2] , self.findLimits()[3] + 0.02 )
             minorLocator    = AutoMinorLocator(4)
             majorLocator    = MultipleLocator(10)
             majorFormatter  = FormatStrFormatter('%f')
             if nrows == ncols == 1:
                 self.axs.grid(linestyle= self.config['grid.linestyle'] , dashes= self.config['grid.dashes'])
                 self.axs.set_title('', color=self.config['title.colors'] )
                 self.axs.xaxis.set_minor_locator(AutoMinorLocator(5))
                 self.axs.yaxis.set_minor_locator(AutoMinorLocator(5))
                 self.axs.yaxis.set_ticks_position('both')
                 self.axs.xaxis.set_ticks_position('both')
                 self.axs.tick_params(which='major', width=self.config['tick.MajW'])  # 1
                 self.axs.tick_params(which='major', length=5.5)
                 self.axs.tick_params(which='minor', width=self.config['tick.minW'])
                 self.axs.tick_params(which='minor', length=3)
                 legend = leg = self.axs.legend()
                 #self.axs.legend( bbox_to_anchor=self.config['legend.bbox'] , loc= self.config['legend.bbox'] )  
                 legend.get_frame().set_linewidth(self.config['legend.framelinewidth'])
                 #legend.get_frame().set_linewidth(0.1)
                 #legend.get_frame().set_ linewidth(0.1)
                 plt.setp(legend.get_texts(), color=self.config['legend.textcolor']  )
             
             elif nrows == 1 or ncols==1 :
                 self.axs[i].grid(linestyle= self.config['grid.linestyle'] , dashes= self.config['grid.dashes'])
                 self.axs[i].set_title('', color='#555555')
                 self.axs[i].xaxis.set_minor_locator(AutoMinorLocator(5))
                 self.axs[i].yaxis.set_minor_locator(AutoMinorLocator(5))
                 self.axs[i].yaxis.set_ticks_position('both')
                 self.axs[i].xaxis.set_ticks_position('both')
                # self.axs[i].tick_params(which='major', width=0.5)
                 self.axs[i].tick_params(which='major', width=self.config['tick.MajW'])  # 1
                 self.axs[i].tick_params(which='major', length=5)
                 self.axs[i].tick_params(which='minor', width=self.config['tick.minW'])
                 self.axs[i].tick_params(which='minor', length=2.5)
                 legend = leg = self.axs[i].legend()
                 legend.get_frame().set_linewidth(0.5)
                 plt.setp(legend.get_texts(), color=self.config['legend.textcolor'])
             else :
                 self.axs[i,j].grid(linestyle= self.config['grid.linestyle'] , dashes= self.config['grid.dashes'])
                 self.axs[i,j].set_title('', color='#555555')
                 self.axs[i,j].xaxis.set_minor_locator(AutoMinorLocator(5))
                 self.axs[i,j].yaxis.set_minor_locator(AutoMinorLocator(5))
                 self.axs[i,j].yaxis.set_ticks_position('both')
                 self.axs[i,j].xaxis.set_ticks_position('both')
                # self.axs[i,j].tick_params(which='major', width=0.5)
                 self.axs[i,j].tick_params(which='major', width=self.config['tick.MajW'])  # 1
                 self.axs[i,j].tick_params(which='major', length=5)
                 self.axs[i,j].tick_params(which='minor', width=self.config['tick.minW'])
                 self.axs[i,j].tick_params(which='minor', length=2.5)
                    
                 #linestyles = ["--","-.",":", (0,(5,2,1,4))]
                # linestyles = [(0,(7,7)),(0,(7,7)),(0,(7,7)), (0,(7,7))]
                # for ax, ls in zip(axs.flat, linestyles):
                #    for spine in ax.spines.values():
                #         spine.set_linestyle(ls)   
                    
                 legend = leg = self.axs[i,j].legend()
             #legend.get_frame().set_linewidth(self.config['legend.framelinewidth'])
                 legend.get_frame().set_linewidth(0.5)
                 plt.setp(legend.get_texts(), color=self.config['legend.textcolor'])
             #axs.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
             ##plt.title(y=1.25)
             #### plt.legend([l1, l2, l3],["HHZ 1", "HHN", "HHE"])
             if nrows == ncols ==1 :
                self.config['subplot.bottom'] = 0.175
                self.config['subplot.right'] = 0.93
                
                plt.subplots_adjust(left=0.125 , bottom=self.config['subplot.bottom'], right=self.config['subplot.right'], top=0.87, wspace=0, hspace=0)
             else:   
                plt.subplots_adjust(left=0.06 , bottom=0.1, right=0.95, top=0.93, wspace=0.2, hspace=0.6) # 0.2 , 0.4
       return self.fig,self.axs

 
    
    
   def __call__ (self,nrows,ncols,figsize=(9.5,4.5)):
       
       self.nrows = nrows 
       self.ncols = ncols 
       self.figsize = figsize
       self.parameters['figure.figsize'] = self.figsize
       self.setparams()
       
       if self.nrows == self.ncols ==1:
            self.fig, self.axs = plt.subplots( self.nrows, self.ncols,figsize=self.figsize)
       else:
            self.fig, self.axs = plt.subplots( self.nrows, self.ncols,figsize=self.figsize)
       
       self.fig,self.axs = self.set( self.nrows, self.ncols,self.fig,self.axs, self.parameters )
       return  self.fig,self.axs

      #  left  = 0.125  # the left side of the subplots of the figure
      #  right = 0.9    # the right side of the subplots of the figure
      #  bottom = 0.1   # the bottom of the subplots of the figure
      #  top = 0.9      # the top of the subplots of the figure
      #  wspace = 0.2   # the amount of width reserved for space between subplots,
      #         # expressed as a fraction of the average axis width
      #  hspace = 0.2   # the amount of height reserved for space between subplots,
      #         # expressed as a fraction of the average axis height
      # 
   


#------------------------------------------------------------------------------------------------

