#! /usr/bin/python
import re
#filename='TCGA-12-0692-01Z-00-DX3_20x_20x_NS-MORPH_1.data'
#markupx=[(8098.367592, 39294.848009)]
#markupy= [(16291.046139,27182.819587)]

#filename= 'TCGA-06-0747-01Z-00-DX2_20x_20x_NS-MORPH_1.data'
#markupx=[(29288.571432,32948.571495)]
#markupy= [(32019.571430,33948.142891)]

filename='TCGA-02-0033-01Z-00-DX1_20x_20x_NS-MORPH_1.data'
markupx=[(26727.655192,30144.896670),(31242.166697,34492.166770),(27544.620709,31168.758744)]
markupy= [(19096.172418,21420.310416),(10488.000009,12896.333397),(2993.034487,5320.620760)]


def compstats():
	f= open(filename,'r');
	i=0;
	xdiff, ydiff=[],[]
	rangex, rangey=(),()

	for line in f:
	#       if(i >1000):
	#               break;
		field= line.rsplit('\"',2)
		pnts= re.findall(r'\d+', field[1])
		x= map(int, pnts[0::2])
		y= map(int, pnts[1::2])
		xdiff.append(max(x)- min(x))
		ydiff.append( max(y)- min(y))
		if i==0:
			rangex=(min(x), max(x))
			rangey=(min(y), max(y))
		else:
			rangex=(min(rangex[0], min(x)), max(rangex[1], max(x)))
			rangey=(min(rangey[0], min(y)), max(rangey[1], max(y)))

		i=i+1
	f.close()
	f= open('stats','w')
	f.write("max_w:{0}\nmax_h:{1}\nmin_w:{2}\nmin_h:{3}\nmax_x:{4}\nmin_x:{5}\nmax_y:{6}\nmin_y:{7}\n".format(max(xdiff), max(ydiff), min(xdiff), min(ydiff), rangex[1],rangex[0], rangey[1], rangey[0]))
	f.close() 

class genimg:
	
	def mark(self, rangex, rangey):
		for x in range(rangex[0]/self.xunit, rangex[1]/self.xunit+1):
			for y in range(rangey[0]/self.yunit, rangey[1]/self.yunit+1):
				self.img[x- self.minxind][y-self.minyind]= self.img[x-self.minxind][y-self.minyind]+1

	def scale(self, newmax):
		m= max([max(x) for x in self.img])
		n= newmax/float(m)
		self.img= [[y*n for y in x] for x in self.img]
		self.img= [map(int, x) for x in self.img]

	def __init__(self, xunit, yunit):
		self.meta={}	
		self.xunit, self.yunit= xunit, yunit
		self.img=[]

		f= open('stats')
		for line in f:
			field= line.split(':')
			self.meta[field[0]]= int(field[1])
		f.close()
		self.minxind= self.meta['min_x']/xunit;
		self.minyind= self.meta['min_y']/yunit;
		
		self.img=[[0 for y in range(self.meta['min_y']/yunit, self.meta['max_y']/yunit+1)] for x in range(self.meta['min_x']/xunit, self.meta['max_x']/xunit+1)]

		f= open(filename,'r');
		i=0;
		for line in f:
		#	if(i >1000):
		#		break;
			field= line.rsplit('\"',2)
			pnts= re.findall(r'\d+', field[1])
			x= map(int, pnts[0::2])
			y= map(int, pnts[1::2])
			self.mark((min(x), max(x)),(min(y), max(y)))
			i=i+1
		self.scale(255)
		f.close()
		for i in range(len(markupx)):
			for x in range(int(markupx[i][0]/xunit), int(markupx[i][1]/xunit)+1):
				self.img[x][int(markupy[i][0]/yunit)]= 255
				self.img[x][int(markupy[i][1]/yunit)]= 255
			for y in range(int(markupy[i][0]/yunit), int(markupy[i][1]/yunit)+1):
				self.img[int(markupx[i][0]/xunit)][y]= 255
		                self.img[int(markupx[i][1]/xunit)][y]= 255
		print '\n'.join([' '.join(map(str, x)) for x in self.img])
compstats()
inst= genimg(100, 80)
