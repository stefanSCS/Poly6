
'''
--This code generates the Poly6-model of a sheet metal based uniaxial and biaxial mechanical test data 
--It is released under the MIT licence by its author Stefan C. Soare
--Read the accompanying 'README.txt' file for instructions
'''

import numpy as np
from matplotlib import pyplot as plt


##global variables and constants
PI=np.pi;PI4=PI/4;PI2=PI/2;PI12=PI/12
oneSix=1.0/6.0
imgFormat='png' ###extension of the image format
missing='*'  ###use this string as placeholder for missing data  
##vCoeff = the vector of Poly6 parameters [1,a_1,a_2,...,a_{16}] 
vCoeff=np.ones(17)
K1=1.0;K2=1.0;Sb=1.0; ##defined in function 'Poly2'
##vTheta = vector of uniaxial directions 
##(modify both material data structure and this vector if other angles are used)
vTheta=[PI12,PI/6,PI/3,5*PI12]
nTheta=len(vTheta)
##vPHI and vTT = vectors of spherical coords on the unit sphere 
vPHI=[PI12/6,PI12/3,2*PI12/3,PI12]
vTT=[-PI/2+k*PI/6 for k in range(1,6)]
nPHI=len(vPHI);nTT=len(vTT)
nP2=nPHI*nTT
##vP2 = the vector of sampling points on Poly2 
vP2=np.zeros((nP2,4))


####-------MATERIAL DATA SAMPLES (use a file instead, or a db, if in a GUI context)-----
##check poly6Params function on isotropic material 
##solution should be:
##a1=1, a2=-3, a3=6, a4=-7, a5=6, a6=-3, a7=1, a8=9, 
##a9=-18, a10=27, a11=-18, a12=9, a13=27, a14=-27, a15=27, a16=27
ISO={
'type':'ISO',
'name':'Mises',
'biax':{'s0':1.0000,'r0':1,'s90':1,'r90':1,'sb':1.0,'rb':1,'bxShape':1.1547},
'UTstress':[1,1,1,1,1],
'UTrvalue':[1,1,1,1,1],
###'PS':[2.0/np.sqrt(3.0),2.0/np.sqrt(3.0)]
}
##-------------------------------------------------------------------------------
AA3104H19={
'type':'AL',
'name':'AA3104-H19',
'biax':{'s0':1.0000,'r0':0.408,'s90':1.051,'r90':1.416,'sb':1.021,'rb':0.359, 'bxShape':'*'},
'UTstress':[ 1.000 , 1.007 , 1.011 , 1.018 , 1.036],
'UTrvalue':[0.475 , 0.639 , 0.984 , 1.060 , 1.173],
'PS':[1.072, 1.239]
}
##-------------------------------------------------------------------------------
AA2024T3={
'type':'AL',
'name':'AA2024-T3',
'biax':{'s0':1.0000,'r0':0.777,'s90':0.954,'r90':0.622,'sb':1.0000,'rb':1.0000, 'bxShape':'*'},
'UTstress':[0.973, 0.939, 0.922, 0.914, 0.929],
'UTrvalue':[0.897, 1.027, 1.075, 0.920, 0.677],
'PS':[1.0753, 1.0493]
}
##--------------------------------------------------------------------------------
AA2090T3={
'type':'AL',
'name':'AA2090-T3',
'biax':{'s0':1.0000,'r0':0.2115,'s90':0.9102,'r90':0.6923,'sb':1.0350,'rb':0.6700,'bxShape':'*'},
'UTstress':[0.9605, 0.9102, 0.8114, 0.8096, 0.8815],
'UTrvalue':[0.3269, 0.6923, 1.5769, 1.0385, 0.5384]
}
##--------------------------------------------------------------------------------
DP980={
'type':'FE',
'name':'DP980',
'biax':{'s0':1.0000,'r0':0.7514,'s90':1.0401, 'r90': 0.7614, 'sb':'*', 'rb':'*','bxShape':'*'},
'UTstress':[ 1.0221 ,  1.0206 ,  1.0103  ,   1.0371 , 1.0352 ],
'UTrvalue':[0.7912 , 0.8784 , 0.8632 , 0.8719 ,	0.8046 ]
}
##--------------------------------------------------------------------------------
AISI304={
'type':'FE',
'name':'AISI-304',
'biax':{'s0':1.0000,'r0':0.822,'s90':0.979, 'r90': 0.798, 'sb':'*', 'rb':'*','bxShape':'*'},
'UTstress':[ '*' ,  '*' ,  0.947  ,   '*' , '*' ],
'UTrvalue':['*' , '*' , 1.104 , '*' ,	'*' ]
}
##--------------------------------------------------------------------------------
#data=AA2090T3
#data=DP980
#data=AA3104H19
#data={}
data=ISO


def rUp(x):
    return float(int(x*10**4+1)/10**4)
    
def rUp6(x):
    return float(int(x*10**6+1)/10**6)    

def rDown(x):
    return float(int(x*10**4)/10**4)    

def dataCheck():
    '''
    --performs basic check on input material data 
    '''
    if(data['biax']):
        if(np.absolute(data['biax']['s0']-1.0)>0.000001):
            print('stress data must be normalized with sigma_0. \nCalculations aborted')
            exit()
    else:
        print('biax-data missing. \nCalculations aborted')
        exit()            
    if(len(data['UTrvalue']) != 5 or len(data['UTstress']) != 5):
        print('incorrect r-value data\nCalculations aborted')
        exit()
    znData=False   
    if(missing in data['UTrvalue'] or missing in data['UTstress']):
        znData=True    
    if(data['biax']['sb']==missing):
        data['biax']['sb']=1.0
    if(data['biax']['rb']==missing):
        data['biax']['rb']=1.0    
    if(data['biax']['bxShape']==missing):
        if(data['type']=='AL'):
            data['biax']['bxShape']=0.01
        if(data['type']=='FE'):
            data['biax']['bxShape']=0.35  
    else:
        if(data['biax']['bxShape']<0.01 or data['biax']['bxShape']>=1.0):
            print('incorrect bxShape provided (must be >0.01 and <1)\nCalculations aborted');exit()        
    zs90=data['biax']['s90'];zs45=data['UTstress'][2];zr45=data['UTrvalue'][2]
    B90min=0.5+(0.5/zs90-1)/zs90; B90max=0.5+(0.5/zs90+1)/zs90; T45=1.0/((1+zr45)*zs45*zs45)
    print('Hill48: B90min = {vv1:.4f}\nHill48: T45 = {vv2:.4f}\nHill48: B90max = {vv3:.4f}'.format(vv1=B90min,vv2=T45,vv3=B90max))
    if(T45<B90min or T45>B90max):
        print('input data s90, s45 and r45 incompatible with Hill48 convexity.\nCalculations aborted')
        exit()
    return [znData,B90min,T45,B90max]



def bxParamPlot():
    '''
    ---calculates the biaxial parameters [a_1,a_2,a_3,a_4,a_5,a_6,a_7] 
    ---plots the biaxial curves correponding to the two extreme values of a_4
    '''
    s0=data['biax']['s0'];r0=data['biax']['r0']
    s90=data['biax']['s90'];r90=data['biax']['r90']
    sb=data['biax']['sb']
    tx=sb/(s0+sb);ty=sb/(s90+sb)
    sb=sb**6;rb=data['biax']['rb']
    vCoeff[2]=(-6.0*r0)/(1.0+r0)
    vCoeff[7]=1.0/(data['biax']['s90'])**6
    vCoeff[6]=(-6.0*r90*vCoeff[7])/(1.0+r90)
    Q1=3.0/((1.0+rb)*sb)+0.5*(vCoeff[6]+2*vCoeff[7]-4*vCoeff[1]-3*vCoeff[2]-2/sb)
    Q2=1.0/sb-(vCoeff[1]+vCoeff[2]+vCoeff[6]+vCoeff[7]+Q1)
    tx2=tx*tx;tx4=tx2*tx2;ty2=ty*ty;ty4=ty2*ty2;tx2ty=tx2+ty4;tx4ty=tx4+ty2
    TT=tx2*tx+ty2*ty-0.5*(tx2ty+tx4ty)
    Q3=vCoeff[1]*(1+ty4*ty2)+vCoeff[2]*(tx+ty4*ty)+Q1*tx2ty+Q2*tx4ty+vCoeff[6]*(tx4*tx+ty)+vCoeff[7]*(1+tx4*tx2)
    domg=PI/100
    vOmega=np.arange(domg,PI4,domg)
    vOmega=np.concatenate((vOmega,np.arange(PI4+domg,PI2,domg)))
    vOmega=np.concatenate((vOmega,np.arange(PI2,PI,domg)))
    #vsMin=np.zeros(vOmega.size);vsMax=np.zeros(vOmega.size)
    c0=0.0;c1=0.0;c2=0.0;sx=0.0;sy=0.0;sx2=0.0;sx4=0.0;sy2=0.0;sy4=0.0;sxsy=0.0;sxsy2=0.0
    s2omg=0.0
    sMin=-10**9;sMax=-sMin;s1=0.0;s2=0.0;dt=0.0;#kk=0
    for omg in vOmega:
        sx=np.cos(omg);sx2=sx*sx;sx4=sx2*sx2
        sy=np.sin(omg);sy2=sy*sy;sy4=sy2*sy2
        sxsy=sx*sy;sxsy2=sxsy*sxsy
        pp=(vCoeff[1]*sx2+vCoeff[2]*sxsy+Q1*sy2)*sx4+(Q2*sx2+vCoeff[6]*sxsy+vCoeff[7]*sy2)*sy4   
        dpdx=(6*vCoeff[1]*sx+5*vCoeff[2]*sy)*sx4+4*Q1*sx*sxsy2+(2*Q2*sx+vCoeff[6]*sy)*sy4
        dpdy=(vCoeff[2]*sx+2*Q1*sy)*sx4+4*Q2*sy*sxsy2+(5*vCoeff[6]*sx+6*vCoeff[7]*sy)*sy4
        dpdxdx=(30*vCoeff[1]*sx2+20*vCoeff[2]*sxsy+12*Q1*sy2)*sx2+2*Q2*sy4
        dpdxdy=(5*vCoeff[2]*sx2+8*Q1*sxsy)*sx2+(8*Q2*sxsy+5*vCoeff[6]*sy2)*sy2
        dpdydy=2*Q1*sx4+(12*Q2*sx2+20*vCoeff[6]*sxsy+30*vCoeff[7]*sy2)*sy2
        dpp=dpdy*sx-dpdx*sy
        ddpp=(dpdxdx*sy-dpdy)*sy+(dpdydy*sx-dpdx)*sx-2*dpdxdy*sxsy
        s2omg=np.sin(2*omg)
        qq=0.125*s2omg*s2omg*(s2omg-1.0)
        dqq=0.125*np.sin(4*omg)*(3*s2omg-2.0)
        ddqq=0.25*(2*np.cos(4*omg)*(3*s2omg-2)+3*np.sin(4*omg)*np.cos(2*omg))        
        #c2=0.125*(4.5*np.sin(2*omg)+3.0*np.cos(4*omg)-1.5*np.sin(6*omg)-0.25*np.cos(8*omg)-2.75)
        c2=6*qq*(6*qq+ddqq)-5*dqq*dqq
        c1=72*pp*qq-10*dpp*dqq+6*(pp*ddqq+ddpp*qq)
        c0=6*pp*(6*pp+ddpp)-5*dpp*dpp
        dt=c1*c1-4*c2*c0
        if(dt<0):
            print('Calculations aborted: data is not consistent with convexity !')
            return             
        dt=np.sqrt(c1*c1-4*c2*c0)
        s1=-0.5*(c1+dt)/c2;s2=0.5*(-c1+dt)/c2
        dt=max(s1,s2);s1=min(s1,s2);s2=dt;
        if(sMin<s1):
            if(s1<sMax):
                sMin=s1
            else:
                print('Calculations aborted: data is not consistent with convexity !')
                return                
        if(sMax>s2):
            if(s2>sMin):
                sMax=s2
            else:
                print('Calculations aborted: data is not consistent with convexity !')
                return                 
        #vsMin[kk]=s1;vsMax[kk]=s2;kk+=1
    A4min=sMin;A4max=sMax    
    print('A4min = ',A4min,'  A4max = ',A4max)   
    dt=(2/(TT*sMax+Q3))**oneSix;sMin=(2/(TT*sMin+Q3))**oneSix;sMax=dt
    print('sMin = {:.4f}'.format(rUp(sMin)),'  sMax = {:.4f}'.format(rDown(sMax)))      
    ##Calculate sxy=constant sections through the yield surface 
    print('Calculating contours...')
    ########generate grid 
    delta = 0.005
    vx = np.arange(-1.25, 1.25, delta)
    vy = np.arange(-1.25, 1.25, delta)
    X, Y = np.meshgrid(vx, vy)
    nGrid=vx.size
    Z=np.zeros((nGrid,nGrid))
    ######prepare the figure for plotting contours
    fig_1=plt.figure()    
    ax_1=fig_1.add_subplot(1,1,1)
    ax_1.set_aspect('equal')
    vCoeff[4]=A4max;vCoeff[3]=Q1-0.5*vCoeff[4];vCoeff[5]=Q2-0.5*vCoeff[4]
    for kk in range(0,nGrid):
        sx=vx[kk];sx2=sx*sx;sx4=sx2*sx2
        for jj in range(0,nGrid):
            sy=vy[jj];sy2=sy*sy;sy4=sy2*sy2;sxsy=sx*sy
            Z[jj,kk]=(vCoeff[1]*sx2+vCoeff[2]*sxsy+vCoeff[3]*sy2)*sx4+vCoeff[4]*sxsy**3+(vCoeff[5]*sx2+vCoeff[6]*sxsy+vCoeff[7]*sy2)*sy4
    CS = ax_1.contour(X, Y, Z,levels=[1.0],colors='k')
    ax_1.grid()
    ax_1.text(x=0.85,y=-1.1,s=r'$\sigma_x/\sigma_0$',fontsize=14)
    ax_1.text(x=-1.15,y=1.0,s=r'$\sigma_y/\sigma_0$',fontsize=14)
    ax_1.text(x=0.02,y=0.04,s=r'$a_4 = A_4^{max}$',fontsize=13)
    fig_1.savefig('{name}_BX_max.{ext}'.format(name=data['name'],ext=imgFormat), bbox_inches='tight')    
    fig_2=plt.figure()    
    ax_2=fig_2.add_subplot(1,1,1)
    ax_2.set_aspect('equal')
    vCoeff[4]=A4min;vCoeff[3]=Q1-0.5*vCoeff[4];vCoeff[5]=Q2-0.5*vCoeff[4]
    for kk in range(0,nGrid):
        sx=vx[kk];sx2=sx*sx;sx4=sx2*sx2
        for jj in range(0,nGrid):
            sy=vy[jj];sy2=sy*sy;sy4=sy2*sy2;sxsy=sx*sy
            Z[jj,kk]=(vCoeff[1]*sx2+vCoeff[2]*sxsy+vCoeff[3]*sy2)*sx4+vCoeff[4]*sxsy**3+(vCoeff[5]*sx2+vCoeff[6]*sxsy+vCoeff[7]*sy2)*sy4
    CS = ax_2.contour(X, Y, Z,levels=[1.0],colors='k')    
    #ax.clabel(CS, inline=True, fontsize=10)
    ax_2.grid()
    ax_2.text(x=0.85,y=-1.1,s=r'$\sigma_x/\sigma_0$',fontsize=14)
    ax_2.text(x=-1.15,y=1.0,s=r'$\sigma_y/\sigma_0$',fontsize=14)  
    ax_2.text(x=0.02,y=0.04,s=r'$a_4 = A_4^{min}$',fontsize=13)
    fig_2.savefig('{name}_BX_min.{ext}'.format(name=data['name'],ext=imgFormat), bbox_inches='tight')
    plt.show()
###-------------END OF 'bxParamPlot' function------------------------------------------------------------------------------------------

def bxParam():
    '''
    ---calculates the biaxial parameters [a_1,a_2,a_3,a_4,a_5,a_6,a_7] 
    '''
    global vCoeff
    s0=data['biax']['s0'];r0=data['biax']['r0']
    s90=data['biax']['s90'];r90=data['biax']['r90']
    sb=data['biax']['sb']
    tx=sb/(s0+sb);ty=sb/(s90+sb)
    sb=sb**6;rb=data['biax']['rb']
    vCoeff[2]=(-6.0*r0)/(1.0+r0)
    vCoeff[7]=1.0/(data['biax']['s90'])**6
    vCoeff[6]=(-6.0*r90*vCoeff[7])/(1.0+r90)
    Q1=3.0/((1.0+rb)*sb)+0.5*(vCoeff[6]+2*vCoeff[7]-4*vCoeff[1]-3*vCoeff[2]-2/sb)
    Q2=1.0/sb-(vCoeff[1]+vCoeff[2]+vCoeff[6]+vCoeff[7]+Q1)
    tx2=tx*tx;tx4=tx2*tx2;ty2=ty*ty;ty4=ty2*ty2;tx2ty=tx2+ty4;tx4ty=tx4+ty2
    TT=tx2*tx+ty2*ty-0.5*(tx2ty+tx4ty)
    Q3=vCoeff[1]*(1+ty4*ty2)+vCoeff[2]*(tx+ty4*ty)+Q1*tx2ty+Q2*tx4ty+vCoeff[6]*(tx4*tx+ty)+vCoeff[7]*(1+tx4*tx2)
    domg=PI/100
    vOmega=np.arange(domg,PI4,domg)
    vOmega=np.concatenate((vOmega,np.arange(PI4+domg,PI2,domg)))
    vOmega=np.concatenate((vOmega,np.arange(PI2,PI,domg)))
    #vsMin=np.zeros(vOmega.size);vsMax=np.zeros(vOmega.size)
    c0=0.0;c1=0.0;c2=0.0;sx=0.0;sy=0.0;sx2=0.0;sx4=0.0;sy2=0.0;sy4=0.0;sxsy=0.0;sxsy2=0.0
    s2omg=0.0
    sMin=-10**9;sMax=-sMin;s1=0.0;s2=0.0;dt=0.0;#kk=0
    for omg in vOmega:
        sx=np.cos(omg);sx2=sx*sx;sx4=sx2*sx2
        sy=np.sin(omg);sy2=sy*sy;sy4=sy2*sy2
        sxsy=sx*sy;sxsy2=sxsy*sxsy
        pp=(vCoeff[1]*sx2+vCoeff[2]*sxsy+Q1*sy2)*sx4+(Q2*sx2+vCoeff[6]*sxsy+vCoeff[7]*sy2)*sy4   
        dpdx=(6*vCoeff[1]*sx+5*vCoeff[2]*sy)*sx4+4*Q1*sx*sxsy2+(2*Q2*sx+vCoeff[6]*sy)*sy4
        dpdy=(vCoeff[2]*sx+2*Q1*sy)*sx4+4*Q2*sy*sxsy2+(5*vCoeff[6]*sx+6*vCoeff[7]*sy)*sy4
        dpdxdx=(30*vCoeff[1]*sx2+20*vCoeff[2]*sxsy+12*Q1*sy2)*sx2+2*Q2*sy4
        dpdxdy=(5*vCoeff[2]*sx2+8*Q1*sxsy)*sx2+(8*Q2*sxsy+5*vCoeff[6]*sy2)*sy2
        dpdydy=2*Q1*sx4+(12*Q2*sx2+20*vCoeff[6]*sxsy+30*vCoeff[7]*sy2)*sy2
        dpp=dpdy*sx-dpdx*sy
        ddpp=(dpdxdx*sy-dpdy)*sy+(dpdydy*sx-dpdx)*sx-2*dpdxdy*sxsy
        s2omg=np.sin(2*omg)
        qq=0.125*s2omg*s2omg*(s2omg-1.0)
        dqq=0.125*np.sin(4*omg)*(3*s2omg-2.0)
        ddqq=0.25*(2*np.cos(4*omg)*(3*s2omg-2)+3*np.sin(4*omg)*np.cos(2*omg))        
        #c2=0.125*(4.5*np.sin(2*omg)+3.0*np.cos(4*omg)-1.5*np.sin(6*omg)-0.25*np.cos(8*omg)-2.75)
        c2=6*qq*(6*qq+ddqq)-5*dqq*dqq
        c1=72*pp*qq-10*dpp*dqq+6*(pp*ddqq+ddpp*qq)
        c0=6*pp*(6*pp+ddpp)-5*dpp*dpp
        dt=c1*c1-4*c2*c0
        if(dt<0):
            print('Calculations aborted: data is not consistent with convexity !')
            return             
        dt=np.sqrt(c1*c1-4*c2*c0)
        s1=-0.5*(c1+dt)/c2;s2=0.5*(-c1+dt)/c2
        dt=max(s1,s2);s1=min(s1,s2);s2=dt;
        if(sMin<s1):
            if(s1<sMax):
                sMin=s1
            else:
                print('Calculations aborted: data is not consistent with convexity !')
                return                
        if(sMax>s2):
            if(s2>sMin):
                sMax=s2
            else:
                print('Calculations aborted: data is not consistent with convexity !')
                return                 
        #vsMin[kk]=s1;vsMax[kk]=s2;kk+=1
    A4min=sMin;A4max=sMax      
    dt=(2/(TT*sMax+Q3))**oneSix;sMin=(2/(TT*sMin+Q3))**oneSix;sMax=dt
    s1=rUp(sMin); s2=rDown(sMax)
    dt=(s1+data['biax']['bxShape']*(s2-s1))**6
    ##if(data['type']=='AL'):
    ##    dt=(s1+0.01*(s2-s1))**6
    ##if(data['type']=='FE'):
    ##    dt=(s1+0.35*(s2-s1))**6
    ##if(data['type']=='ISO'):
    ##    dt=(2/np.sqrt(3.0))**6    
    vCoeff[4]=(2/dt-Q3)/TT;vCoeff[3]=Q1-0.5*vCoeff[4];vCoeff[5]=Q2-0.5*vCoeff[4]
    return((s1,s2,A4min,A4max))    
###-------------END OF 'bxParam' function------------------------------------------------------------------------------------------    

def cvxCheck():
    '''
    performs convexity check of Poly6 for a set of coefficients
    returns 'True' if convex, 'False' otherwise
    '''
    vvGrid=100
    dTheta=(2*PI)/(4*vvGrid); dPhi=PI/(2*vvGrid)
    h11=0.0;h12=0.0;h13=0.0;h22=0.0;h23=0.0;h33=0.0;det2=0.0;det3=0.0;
    dPdx=0.0;dPdy=0.0;dPdxy=0.0;PP=0.0;fiveSix=5.0/6.0
    if(4*vCoeff[13]*vCoeff[15]-vCoeff[14]*vCoeff[14]<0): ##the case PHI=0
        print('not convex at top')
        return(False)
    PHI=dPhi    
    convex=True; exitFlag=False 
    while(PHI<PI2):
        sxy=np.cos(PHI);sxy2=sxy*sxy;sxy3=sxy*sxy2;sxy4=sxy*sxy3;sxy6=sxy4*sxy2;sp=np.sin(PHI);PHI+=dPhi
        THETA=0.0
        while(THETA<PI):
            sx=sp*np.cos(THETA);sy=sp*np.sin(THETA);THETA+=dTheta
            sx2=sx*sx;sx3=sx*sx2;sx4=sx3*sx;sxsy=sx*sy;sxsy2=sxsy*sxsy;sy2=sy*sy;sy3=sy*sy2;sy4=sy3*sy
            PP=(vCoeff[1]*sx2+vCoeff[2]*sxsy+vCoeff[3]*sy2)*sx4+vCoeff[4]*sxsy**3+\
            (vCoeff[5]*sx2+vCoeff[6]*sxsy+vCoeff[7]*sy2)*sy4+\
            (vCoeff[8]*sx4+(vCoeff[9]*sx2+vCoeff[10]*sxsy+vCoeff[11]*sy2)*sxsy+vCoeff[12]*sy4)*sxy2+\
            (vCoeff[13]*sx2+vCoeff[14]*sxsy+vCoeff[15]*sy2)*sxy4+vCoeff[16]*sxy6  
            if(rUp6(PP)<0):
                convex=False;exitFlag=True; print('PP negative: ', PP)
                break            
            dPdx=(6*vCoeff[1]*sx2+5*vCoeff[2]*sxsy+4*vCoeff[3]*sy2)*sx3+\
            (3*vCoeff[4]*sx2+2*vCoeff[5]*sxsy+vCoeff[6]*sy2)*sy3+\
            (4*vCoeff[8]*sx3+(3*vCoeff[9]*sx+2*vCoeff[10]*sy)*sxsy+vCoeff[11]*sy3)*sxy2+\
            (2*vCoeff[13]*sx+vCoeff[14]*sy)*sxy4
            h11=10*(3*vCoeff[1]*sx2+2*vCoeff[2]*sxsy)*sx2+12*vCoeff[3]*sxsy2+(6*vCoeff[4]*sxsy+2*vCoeff[5]*sy2)*sy2+\
            (12*vCoeff[8]*sx2+6*vCoeff[9]*sxsy+2*vCoeff[10]*sy2)*sxy2+2*vCoeff[13]*sxy4
            h11-=fiveSix*dPdx*dPdx/PP
            if(rUp6(h11)<0):
                convex=False;exitFlag=True; 
                print('h11 negative: ', h11)
                print(PHI,THETA);#exit()
                break            
            dPdy=(vCoeff[2]*sx2+2*vCoeff[3]*sxsy+3*vCoeff[4]*sy2)*sx3+\
            (4*vCoeff[5]*sx2+5*vCoeff[6]*sxsy+6*vCoeff[7]*sy2)*sy3+\
            (vCoeff[9]*sx3+(2*vCoeff[10]*sx+3*vCoeff[11]*sy)*sxsy+4*vCoeff[12]*sy3)*sxy2+\
            (vCoeff[14]*sx+2*vCoeff[15]*sy)*sxy4
            dPdxy=2*((vCoeff[8]*sx+vCoeff[9]*sy)*sx3+vCoeff[10]*sx2*sy2+(vCoeff[11]*sx+vCoeff[12]*sy)*sy3)*sxy+\
            4*(vCoeff[13]*sx2+vCoeff[14]*sxsy+vCoeff[15]*sy2)*sxy3+6*vCoeff[16]*sxy3*sxy2        
            h12=(5*vCoeff[2]*sx2+8*vCoeff[3]*sxsy)*sx2+9*vCoeff[4]*sxsy2+(8*vCoeff[5]*sxsy+5*vCoeff[6]*sy2)*sy2+\
            (3*vCoeff[9]*sx2+4*vCoeff[10]*sxsy+3*vCoeff[11]*sy2)*sxy2+vCoeff[14]*sxy4
            h12-=fiveSix*dPdx*dPdy/PP
            h22=(2*vCoeff[3]*sx2+6*vCoeff[4]*sxsy)*sx2+12*vCoeff[5]*sxsy2+10*(2*vCoeff[6]*sxsy+3*vCoeff[7]*sy2)*sy2+\
            2*(vCoeff[10]*sx2+3*vCoeff[11]*sxsy+6*vCoeff[12]*sy2)*sxy2+2*vCoeff[15]*sxy4
            h22-=fiveSix*dPdy*dPdy/PP
            det2=h11*h22-h12*h12
            if(rUp6(det2)<0):
                convex=False;exitFlag=True; print('det2 negative: ', det2)
                break
            h13=2*((4*vCoeff[8]*sx+3*vCoeff[9]*sy)*sx2+(2*vCoeff[10]*sx+vCoeff[11]*sy)*sy2)*sxy+\
            4*(2*vCoeff[13]*sx+vCoeff[14]*sy)*sxy3
            h13-=fiveSix*dPdx*dPdxy/PP
            h23=2*((vCoeff[9]*sx+2*vCoeff[10]*sy)*sx2+(3*vCoeff[11]*sx+4*vCoeff[12]*sy)*sy2)*sxy+\
            4*(vCoeff[14]*sx+2*vCoeff[15]*sy)*sxy3
            h23-=fiveSix*dPdy*dPdxy/PP
            h33=2*((vCoeff[8]*sx2+vCoeff[9]*sxsy)*sx2+vCoeff[10]*sxsy2+(vCoeff[11]*sxsy+vCoeff[12]*sy2)*sy2)+\
            12*(vCoeff[13]*sx2+vCoeff[14]*sxsy+vCoeff[15]*sy2)*sxy2+30*vCoeff[16]*sxy4
            h33-=fiveSix*dPdxy*dPdxy/PP
            det3=h33*det2-h23*(h11*h23-h12*h13)+h13*(h12*h23-h13*h22)
            if(rUp6(det3)<0):
                convex=False;exitFlag=True; print('det3 negative: ', det3)
                break
        if(exitFlag):
            break
    return convex
###-------------END OF 'cvxCheck' function------------------------------------------------------------------------------------------


def Poly2():
    '''
    calculates the Poly2 coeffs, a_{16} and the Poly2 sampling points
    '''
    global K1, K2, Sb
    s00=data['biax']['s0'];s90=data['biax']['s90']
    s45=data['UTstress'][2]
    r45=data['UTrvalue'][2]    
    ss45=2*s00/s45;ss45=ss45*ss45;ss45=ss45**3
    Sb=vCoeff[1]+vCoeff[2]+vCoeff[3]+vCoeff[4]+vCoeff[5]+vCoeff[6]+vCoeff[7]
    K1=2*Sb+ss45*(r45-0.5)/(r45+1)
    K2=ss45*(2*r45+0.5)/(r45+1)+Sb
    HA1=s00;HA3=(s00/s90)**2;HA4=((r45+0.5)/(r45+1))*(2*s00/s45)**2
    HA2=(2*s00/s45)**2-(HA1+HA3+HA4)
    vCoeff[16]=HA4**3
    ii=0;rho=0.0
    for kk in range(0,nPHI):
        PHI=vPHI[kk];cp=np.cos(PHI);sp=np.sin(PHI);cp2=cp*cp;sp2=sp*sp
        for jj in range(0,nTT):
            TT=vTT[jj];ct=np.cos(TT);st=np.sin(TT)
            rho=1/np.sqrt((HA1*ct*ct+HA2*st*ct+HA3*st*st)*sp2+HA4*cp2)
            vP2[ii]=[sp*ct,sp*st,cp,rho]
            ii+=1
###-------------END OF 'Poly2' function------------------------------------------------------------------------------------------

def uaxParam():
    '''
    calculates the Poly6 parameters a_8,..., a_{15} 
    '''
    global vCoeff
    nEq=2*nTheta+nP2; #nEq=nTheta
    matA=np.zeros((nEq,6))
    vB=np.zeros(nEq)
    vStress=np.ones(4)
    vStress[0]=data['UTstress'][0];vStress[1]=data['UTstress'][1]
    vStress[2]=data['UTstress'][3];vStress[3]=data['UTstress'][4]
    vRval=np.zeros(4)
    vRval[0]=data['UTrvalue'][0];vRval[1]=data['UTrvalue'][1]
    vRval[2]=data['UTrvalue'][3];vRval[3]=data['UTrvalue'][4] 
    jj=0;BB=0.0;CC=0.0;DD=0.0
    for kk in range(0,nTheta):  ## calculate the matrix entries from uniaxial data 
        theta=vTheta[kk]; cc=np.cos(theta);ss=np.sin(theta)
        x=cc*cc;y=ss*ss;z=cc*ss
        x2=x*x;x3=x2*x;x4=x3*x;y2=y*y;y3=y2*y;y4=y3*y;xy=x*y;xmy=x-y;z2=z*z;z4=z2*z2
        matA[kk][0]=x2
        matA[kk][1]=x2*y
        matA[kk][2]=-x*y2
        matA[kk][3]=-y2
        matA[kk][4]=x*z2
        matA[kk][5]=-y*z2 
        BB=(vCoeff[1]*x2+vCoeff[2]*xy+vCoeff[3]*y2)*x2*x2+vCoeff[4]*xy*xy*xy+(vCoeff[5]*x2+vCoeff[6]*xy+vCoeff[7]*y2)*y2*y2 
        st=(1.0/vStress[kk])**6
        vB[kk]=((st-BB)/z2+z4*(K1-K2))/xmy
        ##continue
        jj=nTheta+kk
        rr=vRval[kk]
        matA[jj][0]=2*((2*x3-z2)*rr-x4+(2*x2-1+3*z2)*z2)  #; matA[jj][0]*=z2
        matA[jj][1]=(x3+(3*x-2)*z2)*rr+x4+(-2*x2-2+9*z2)*z2  #; matA[jj][1]*=z2
        matA[jj][2]=(y3+(3*y-2)*z2)*rr+y4+(-2*y2-2+9*z2)*z2  #; matA[jj][2]*=z2
        matA[jj][3]=2*((2*y3-z2)*rr-y4+(2*y2-1+3*z2)*z2)  #; matA[jj][3]*=z2
        matA[jj][4]=z2*((2*x-1)*rr+8*z2-4*x2-1)   #; matA[jj][4]*=z2
        matA[jj][5]=z2*((2*y-1)*rr+8*z2-4*y2-1)   #; matA[jj][5]*=z2
        CC=(rr+y)*((6*vCoeff[1]*x2+5*vCoeff[2]*xy+4*vCoeff[3]*y2)*x3+(3*vCoeff[4]*x2+2*vCoeff[5]*xy+vCoeff[6]*y2)*y3)
        CC+=(rr+x)*((vCoeff[2]*x2+2*vCoeff[3]*xy+3*vCoeff[4]*y2)*x3+(4*vCoeff[5]*x2+5*vCoeff[6]*xy+6*vCoeff[7]*y2)*y3)
        DD=K2*(rr+1-6*z2)-2*K1*(rr+1-3*z2)
        vB[jj]=-(DD*z2+CC/z2)
    for jj in range(0,nP2): ## calculate the matrix entries from Poly2 sampling data 
        x,y,z,rho=vP2[jj]
        x2=x*x;y2=y*y;z2=z*z;xy=x*y
        kk=2*nTheta+jj
        matA[kk][0]=x2*(x+y)
        matA[kk][1]=x2*y
        matA[kk][2]=-x*y2 
        matA[kk][3]=-y2*(x+y)
        matA[kk][4]=x*z2
        matA[kk][5]=-y*z2
        BB=(vCoeff[1]*x2+vCoeff[2]*xy+vCoeff[3]*y2)*x2*x2+vCoeff[4]*xy*xy*xy+(vCoeff[5]*x2+vCoeff[6]*xy+vCoeff[7]*y2)*y2*y2 
        vB[kk]=((1/rho**6-BB)/z2+xy*(xy*K1-z2*K2)-vCoeff[16]*(xy-z2)**2)/(x-y)
    MAT=np.zeros((6,6),dtype=np.float64) ##final system matrix 
    YY=np.zeros(6,dtype=np.float64) ##right-hand side 
    vDiag=np.ones(nEq,dtype=np.float64) ## the vector of weights 
    convex=False
    LL1=1.0000;LL2=LL1;dL=0.1;  ss=0; tt=0
    while(LL2>-0.0001): ## find first L for which Poly6 is convex 
        vDiag[0:2*nTheta]=LL2/(2*nTheta); vDiag[2*nTheta:]=(1.0-LL2)/nP2
        for ii in range(0,6):
            for jj in range(0,6):
                ss=0
                for kk in range(0,nEq):
                    ss+=matA[kk][ii]*matA[kk][jj]*vDiag[kk]
                MAT[ii][jj]=ss  
            tt=0
            for kk in range(0,nEq):
                tt+=matA[kk][ii]*vB[kk]*vDiag[kk]
            YY[ii]=tt                
        vX=np.linalg.solve(MAT,YY)
        vCoeff[8]=vX[0]; vCoeff[9]=vX[1]
        vCoeff[10]=vCoeff[16]-(K1+vX[0]+vX[1]+vX[2]+vX[3])
        vCoeff[11]=vX[2]; vCoeff[12]=vX[3]; vCoeff[13]=vX[4]
        vCoeff[14]=K2-(2*vCoeff[16]+vX[4]+vX[5])
        vCoeff[15]=vX[5] ##print(vCoeff)
        if(cvxCheck()):
            convex=True; break
        LL1=LL2; LL2-=dL
    NLL=0        
    if(convex): ##found an LL for which Poly6 is convex 
        LL=0.0
        while(LL1-LL2>0.0001):##find the best LL via bisection method
            LL=0.5*(LL1+LL2); NLL+=1
            vDiag[0:2*nTheta]=LL/(2*nTheta); vDiag[2*nTheta:]=(1.0-LL)/nP2
            for ii in range(0,6):
                for jj in range(0,6):
                    ss=0
                    for kk in range(0,nEq):
                        ss+=matA[kk][ii]*matA[kk][jj]*vDiag[kk]
                    MAT[ii][jj]=ss  
                tt=0
                for kk in range(0,nEq):
                    tt+=matA[kk][ii]*vB[kk]*vDiag[kk]
                YY[ii]=tt                
            vX=np.linalg.solve(MAT,YY)
            vCoeff[8]=vX[0]; vCoeff[9]=vX[1]
            vCoeff[10]=vCoeff[16]-(K1+vX[0]+vX[1]+vX[2]+vX[3])
            vCoeff[11]=vX[2]; vCoeff[12]=vX[3]; vCoeff[13]=vX[4]
            vCoeff[14]=K2-(2*vCoeff[16]+vX[4]+vX[5])
            vCoeff[15]=vX[5]
            if(cvxCheck()):
                LL2=LL
            else:
                LL1=LL
    return (convex,rDown(LL2),LL1,NLL)         
###-------------END OF 'uaxParam' function------------------------------------------------------------------------------------------

def uaxParamMissing():
    '''
    calculates the Poly6 parameters a_8,..., a_{15} when the directional data set is minimal 
    '''
    global vCoeff
    nEq=nP2
    matA=np.zeros((nEq,6))
    vB=np.zeros(nEq)
    for jj in range(0,nP2): ## calculate the matrix entries from Poly2 sampling data 
        x,y,z,rho=vP2[jj]
        x2=x*x;y2=y*y;z2=z*z;xy=x*y
        kk=2*nTheta+jj
        matA[jj][0]=x2*(x+y)
        matA[jj][1]=x2*y
        matA[jj][2]=-x*y2 
        matA[jj][3]=-y2*(x+y)
        matA[jj][4]=x*z2
        matA[jj][5]=-y*z2
        BB=(vCoeff[1]*x2+vCoeff[2]*xy+vCoeff[3]*y2)*x2*x2+vCoeff[4]*xy*xy*xy+(vCoeff[5]*x2+vCoeff[6]*xy+vCoeff[7]*y2)*y2*y2 
        vB[jj]=((1/rho**6-BB)/z2+xy*(xy*K1-z2*K2)-vCoeff[16]*(xy-z2)**2)/(x-y)
    MAT=np.zeros((6,6),dtype=np.float64) ##final system matrix 
    YY=np.zeros(6,dtype=np.float64) ##right-hand side 
    vDiag=np.ones(nEq,dtype=np.float64) ## the vector of weights 
    convex=False
    vDiag[0:]=1.0/nP2
    for ii in range(0,6):
        for jj in range(0,6):
            ss=0
            for kk in range(0,nEq):
                ss+=matA[kk][ii]*matA[kk][jj]*vDiag[kk]
            MAT[ii][jj]=ss  
        tt=0
        for kk in range(0,nEq):
            tt+=matA[kk][ii]*vB[kk]*vDiag[kk]
        YY[ii]=tt                
    vX=np.linalg.solve(MAT,YY)
    vCoeff[8]=vX[0]; vCoeff[9]=vX[1]
    vCoeff[10]=vCoeff[16]-(K1+vX[0]+vX[1]+vX[2]+vX[3])
    vCoeff[11]=vX[2]; vCoeff[12]=vX[3]; vCoeff[13]=vX[4]
    vCoeff[14]=K2-(2*vCoeff[16]+vX[4]+vX[5])
    vCoeff[15]=vX[5] ##print(vCoeff)
    if(cvxCheck()):
        convex=True
    return (convex,0.0,0.0,0)         
###-------------END OF 'uaxParamMissing' function------------------------------------------------------------------------------------------


def uaxParamAddedData(addedData):
    '''
    calculates the Poly6 parameters a_8,..., a_{15} 
    when data is missing: 'euristic' data is provided instead 
    '''
    global vCoeff
    nnTheta=len(addedData)
    nEq=2*nnTheta+nP2
    matA=np.zeros((nEq,6))
    vB=np.zeros(nEq)
    jj=0;BB=0.0;CC=0.0;DD=0.0
    for kk in range(0,nnTheta):  ## calculate the matrix entries from uniaxial data 
        theta=addedData[kk][0]; cc=np.cos(theta);ss=np.sin(theta)
        x=cc*cc;y=ss*ss;z=cc*ss
        x2=x*x;x3=x2*x;x4=x3*x;y2=y*y;y3=y2*y;y4=y3*y;xy=x*y;xmy=x-y;z2=z*z;z4=z2*z2
        matA[kk][0]=x2
        matA[kk][1]=x2*y
        matA[kk][2]=-x*y2
        matA[kk][3]=-y2
        matA[kk][4]=x*z2
        matA[kk][5]=-y*z2 
        BB=(vCoeff[1]*x2+vCoeff[2]*xy+vCoeff[3]*y2)*x2*x2+vCoeff[4]*xy*xy*xy+(vCoeff[5]*x2+vCoeff[6]*xy+vCoeff[7]*y2)*y2*y2 
        st=(1.0/addedData[kk][1])**6
        vB[kk]=((st-BB)/z2+z4*(K1-K2))/xmy
        ##continue
        jj=nnTheta+kk
        rr=addedData[kk][2]
        matA[jj][0]=2*((2*x3-z2)*rr-x4+(2*x2-1+3*z2)*z2)  #; matA[jj][0]*=z2
        matA[jj][1]=(x3+(3*x-2)*z2)*rr+x4+(-2*x2-2+9*z2)*z2  #; matA[jj][1]*=z2
        matA[jj][2]=(y3+(3*y-2)*z2)*rr+y4+(-2*y2-2+9*z2)*z2  #; matA[jj][2]*=z2
        matA[jj][3]=2*((2*y3-z2)*rr-y4+(2*y2-1+3*z2)*z2)  #; matA[jj][3]*=z2
        matA[jj][4]=z2*((2*x-1)*rr+8*z2-4*x2-1)   #; matA[jj][4]*=z2
        matA[jj][5]=z2*((2*y-1)*rr+8*z2-4*y2-1)   #; matA[jj][5]*=z2
        CC=(rr+y)*((6*vCoeff[1]*x2+5*vCoeff[2]*xy+4*vCoeff[3]*y2)*x3+(3*vCoeff[4]*x2+2*vCoeff[5]*xy+vCoeff[6]*y2)*y3)
        CC+=(rr+x)*((vCoeff[2]*x2+2*vCoeff[3]*xy+3*vCoeff[4]*y2)*x3+(4*vCoeff[5]*x2+5*vCoeff[6]*xy+6*vCoeff[7]*y2)*y3)
        DD=K2*(rr+1-6*z2)-2*K1*(rr+1-3*z2)
        vB[jj]=-(DD*z2+CC/z2)
    for jj in range(0,nP2): ## calculate the matrix entries from Poly2 sampling data 
        x,y,z,rho=vP2[jj]
        x2=x*x;y2=y*y;z2=z*z;xy=x*y
        kk=2*nnTheta+jj
        matA[kk][0]=x2*(x+y)
        matA[kk][1]=x2*y
        matA[kk][2]=-x*y2 
        matA[kk][3]=-y2*(x+y)
        matA[kk][4]=x*z2
        matA[kk][5]=-y*z2
        BB=(vCoeff[1]*x2+vCoeff[2]*xy+vCoeff[3]*y2)*x2*x2+vCoeff[4]*xy*xy*xy+(vCoeff[5]*x2+vCoeff[6]*xy+vCoeff[7]*y2)*y2*y2 
        vB[kk]=((1/rho**6-BB)/z2+xy*(xy*K1-z2*K2)-vCoeff[16]*(xy-z2)**2)/(x-y)
    MAT=np.zeros((6,6),dtype=np.float64) ##final system matrix 
    YY=np.zeros(6,dtype=np.float64) ##right-hand side 
    vDiag=np.ones(nEq,dtype=np.float64) ## the vector of weights 
    convex=False
    LL1=1.0000;LL2=LL1;dL=0.1;  ss=0; tt=0
    while(LL2>-0.0001): ## find first L for which Poly6 is convex 
        vDiag[0:2*nnTheta]=LL2/(2*nnTheta); vDiag[2*nnTheta:]=(1.0-LL2)/nP2
        for ii in range(0,6):
            for jj in range(0,6):
                ss=0
                for kk in range(0,nEq):
                    ss+=matA[kk][ii]*matA[kk][jj]*vDiag[kk]
                MAT[ii][jj]=ss  
            tt=0
            for kk in range(0,nEq):
                tt+=matA[kk][ii]*vB[kk]*vDiag[kk]
            YY[ii]=tt                
        vX=np.linalg.solve(MAT,YY)
        vCoeff[8]=vX[0]; vCoeff[9]=vX[1]
        vCoeff[10]=vCoeff[16]-(K1+vX[0]+vX[1]+vX[2]+vX[3])
        vCoeff[11]=vX[2]; vCoeff[12]=vX[3]; vCoeff[13]=vX[4]
        vCoeff[14]=K2-(2*vCoeff[16]+vX[4]+vX[5])
        vCoeff[15]=vX[5] ##print(vCoeff)
        if(cvxCheck()):
            convex=True; break
        LL1=LL2; LL2-=dL
    NLL=0        
    if(convex): ##found an LL for which Poly6 is convex 
        LL=0.0
        while(LL1-LL2>0.0001):##find the best LL via bisection method
            LL=0.5*(LL1+LL2); NLL+=1
            vDiag[0:2*nnTheta]=LL/(2*nnTheta); vDiag[2*nnTheta:]=(1.0-LL)/nP2
            for ii in range(0,6):
                for jj in range(0,6):
                    ss=0
                    for kk in range(0,nEq):
                        ss+=matA[kk][ii]*matA[kk][jj]*vDiag[kk]
                    MAT[ii][jj]=ss  
                tt=0
                for kk in range(0,nEq):
                    tt+=matA[kk][ii]*vB[kk]*vDiag[kk]
                YY[ii]=tt                
            vX=np.linalg.solve(MAT,YY)
            vCoeff[8]=vX[0]; vCoeff[9]=vX[1]
            vCoeff[10]=vCoeff[16]-(K1+vX[0]+vX[1]+vX[2]+vX[3])
            vCoeff[11]=vX[2]; vCoeff[12]=vX[3]; vCoeff[13]=vX[4]
            vCoeff[14]=K2-(2*vCoeff[16]+vX[4]+vX[5])
            vCoeff[15]=vX[5]
            if(cvxCheck()):
                LL2=LL
            else:
                LL1=LL
    return (convex,rDown(LL2),LL1,NLL)         
###-------------END OF 'uaxParamAddedData' function------------------------------------------------------------------------------------------

def plotPoly6(missingData,addedData,strOption):
    '''
    plots sigma_{xy}=const sections and the directional Poly6-predictions 
    '''
    ##Calculate sxy=constant sections through the yield surface 
    print('Calculating contours...')
    ########generate grid 
    delta = 0.01
    vx = np.arange(-1.25, 1.25, delta)
    vy = np.arange(-1.25, 1.25, delta)
    X, Y = np.meshgrid(vx, vy)
    nGrid=vx.size
    Z=np.zeros((nGrid,nGrid))
    ######prepare the figure for plotting contours
    fig_1=plt.figure()    
    ax_1=fig_1.add_subplot(1,1,1)
    ax_1.set_aspect('equal')
    ######compute level sets and add them to figure 
    maxSXY=1.0/vCoeff[16]**(1/6); deltaSXY=maxSXY/20
    vLevels=[k*deltaSXY for k in range(0,20)]+[maxSXY-0.125*deltaSXY]
    kLevel=1
    for sxy in vLevels:
        sxy2=sxy*sxy;sxy4=sxy2*sxy2;sxy6=sxy4*sxy2
        for kk in range(0,nGrid):
            sx=vx[kk];sx2=sx*sx;sx4=sx2*sx2
            for jj in range(0,nGrid):
                sy=vy[jj];sy2=sy*sy;sy4=sy2*sy2;sxsy=sx*sy
                Z[jj,kk]=(vCoeff[1]*sx2+vCoeff[2]*sxsy+vCoeff[3]*sy2)*sx4+vCoeff[4]*sxsy**3+\
                (vCoeff[5]*sx2+vCoeff[6]*sxsy+vCoeff[7]*sy2)*sy4+\
                (vCoeff[8]*sx4+(vCoeff[9]*sx2+vCoeff[10]*sxsy+vCoeff[11]*sy2)*sxsy+vCoeff[12]*sy4)*sxy2+\
                (vCoeff[13]*sx2+vCoeff[14]*sxsy+vCoeff[15]*sy2)*sxy4+vCoeff[16]*sxy6
        CS = ax_1.contour(X, Y, Z,levels=[1.0],colors='k')
        #ax.clabel(CS, inline=True, fontsize=10)
        print('Plotting level curve {level}/21: sxy = {shear} ...Done'.format(level=str(kLevel),shear="{:.3f}".format(sxy)))
        kLevel+=1
    ax_1.grid()
    ax_1.text(x=0.875,y=-1.1,s=r'$\sigma_x/\sigma_0$',fontsize=14)
    ax_1.text(x=-1.15,y=1.0,s=r'$\sigma_y/\sigma_0$',fontsize=14)
    fig_1.savefig('{name}_contours{opt}.{ext}'.format(name=data['name'],ext=imgFormat,opt=strOption), bbox_inches='tight')
    print('--------------------')
    ##Calculate predicted uniaxial traction data 
    print('Plotting uniaxial data...')
    #######prepare the figure for plotting uniaxial data
    fig_2=plt.figure()    
    ax_21=fig_2.add_subplot(2,1,1)
    ax_22=fig_2.add_subplot(2,1,2)
    nPoints=100
    dTheta=(PI/2)/nPoints;radian=180/PI
    theta=0;kk=0
    vTheta=np.zeros(nPoints+1)
    vStheta=np.zeros(nPoints+1)
    vRtheta=np.zeros(nPoints+1)
    while(kk<nPoints+1):
        vTheta[kk]=theta*radian
        cc=np.cos(theta);ss=np.sin(theta);sx=cc*cc;sy=ss*ss;sxy=cc*ss
        sx2=sx*sx;sx4=sx2*sx2;sy2=sy*sy;sy4=sy2*sy2;sxsy=sx*sy
        sxy2=sxy*sxy;sxy4=sxy2*sxy2;sxy6=sxy4*sxy2
        sTheta=(vCoeff[1]*sx2+vCoeff[2]*sxsy+vCoeff[3]*sy2)*sx4+vCoeff[4]*sxsy**3+\
        (vCoeff[5]*sx2+vCoeff[6]*sxsy+vCoeff[7]*sy2)*sy4+\
        (vCoeff[8]*sx4+(vCoeff[9]*sx2+vCoeff[10]*sxsy+vCoeff[11]*sy2)*sxsy+vCoeff[12]*sy4)*sxy2+\
        (vCoeff[13]*sx2+vCoeff[14]*sxsy+vCoeff[15]*sy2)*sxy4+vCoeff[16]*sxy6
        vStheta[kk]=1/sTheta**oneSix
        sx3=sx2*sx;sy3=sy2*sy;sxy3=sxy2*sxy
        dPdx=(6*vCoeff[1]*sx2+5*vCoeff[2]*sxsy+4*vCoeff[3]*sy2)*sx3+\
        (3*vCoeff[4]*sx2+2*vCoeff[5]*sxsy+vCoeff[6]*sy2)*sy3+\
        (4*vCoeff[8]*sx3+(3*vCoeff[9]*sx+2*vCoeff[10]*sy)*sxsy+vCoeff[11]*sy3)*sxy2+\
        (2*vCoeff[13]*sx+vCoeff[14]*sy)*sxy4
        dPdy=(vCoeff[2]*sx2+2*vCoeff[3]*sxsy+3*vCoeff[4]*sy2)*sx3+\
        (4*vCoeff[5]*sx2+5*vCoeff[6]*sxsy+6*vCoeff[7]*sy2)*sy3+\
        (vCoeff[9]*sx3+(2*vCoeff[10]*sx+3*vCoeff[11]*sy)*sxsy+4*vCoeff[12]*sy3)*sxy2+\
        (vCoeff[14]*sx+2*vCoeff[15]*sy)*sxy4
        dPdxy=2*((vCoeff[8]*sx+vCoeff[9]*sy)*sx3+vCoeff[10]*sx2*sy2+(vCoeff[11]*sx+vCoeff[12]*sy)*sy3)*sxy+\
        4*(vCoeff[13]*sx2+vCoeff[14]*sxsy+vCoeff[15]*sy2)*sxy3+6*vCoeff[16]*sxy3*sxy2
        vRtheta[kk]=(sx*dPdx+sxy*dPdxy+sy*dPdy)/(dPdx+dPdy)-1.0
        theta+=dTheta; kk+=1
    ax_21.plot(vTheta,vStheta,color='k',linewidth=1,label='Prediction')
    ax_21.grid()
    vAddedDataX=[]; vAddedDataY=[]; nAddedData=len(addedData)
    if(missingData):
        vDataX=[0,45,90]
        vDataY=[data['biax']['s0'],data['UTstress'][2],data['biax']['s90']]
        vAddedDataX=[radian*addedData[kk][0] for kk in range(0,nAddedData)]
        vAddedDataY=[addedData[kk][1] for kk in range(0,nAddedData)]
    else:    
        vDataX=[0,15,30,45,60,75,90]
        vDataY=[data['biax']['s0']]+data['UTstress']+[data['biax']['s90']]
    ax_21.plot(vDataX,vDataY,'bo',markerfacecolor='b',markersize=5,label='Exp data')
    ax_21.plot(vAddedDataX,vAddedDataY,'ko',markerfacecolor='w',markersize=6)
    ax_21.set_xticks([0,15,30,45,60,75,90])
    ax_21.set_xticklabels(['0','15','30','45','60','75','90'],fontsize=9)
    maxS=max(vStheta);minS=min(vStheta)
    ax_21.set_ylim([0.975*minS,1.05*maxS])
    ax_21.text(x=0,y=1.01*maxS,s=r'$\sigma_{\theta}$',fontsize=14)
    ax_21.legend(loc=0)
    ax_22.plot(vTheta,vRtheta,color='k',linewidth=1,label='Prediction')
    ax_22.grid()
    if(missingData):
        vDataY=[data['biax']['r0'],data['UTrvalue'][2],data['biax']['r90']]
        vAddedDataY=[addedData[kk][2] for kk in range(0,nAddedData)]
    else:    
        vDataY=[data['biax']['r0']]+data['UTrvalue']+[data['biax']['r90']]
    ax_22.plot(vDataX,vDataY,'bs',markerfacecolor='b',markersize=5,label='Exp data')
    ax_22.plot(vAddedDataX,vAddedDataY,'ko',markerfacecolor='w',markersize=6)
    ax_22.set_xticks([0,15,30,45,60,75,90])
    ax_22.set_xticklabels(['0','15','30','45','60','75','90'],fontsize=9)
    maxS=max(vRtheta);minS=min(vRtheta)
    ax_22.set_ylim([0.87*minS,1.085*maxS])
    ax_22.text(x=0,y=vDataY[0]+0.2*(maxS-vDataY[0]),s=r'$r_{\theta}$',fontsize=14)
    ax_22.text(x=91.5,y=0.95*minS,s=r'$\theta$',fontsize=14)
    ax_22.legend(loc=0)
    fig_2.savefig('{name}_Uniaxial{opt}.{ext}'.format(name=data['name'],ext=imgFormat,opt=strOption), bbox_inches='tight')
    print('--------------------')
    print('Done')
    plt.show()
###-------------END OF 'plotPoly6' function------------------------------------------------------------------------------------------

def readData(fName):
    try:
        ff=open(fName,'r')
    except IOError as err:
        print(err)
        exit()        
    vLine=[];fVal=0.0
    vCheck={'name':False,'type':False,'option':False,
            's0':False,'s15':False,'s30':False,'s45':False,'s60':False,'s75':False,'s90':False,
            'r0':False,'r15':False,'r30':False,'r45':False,'r60':False,'r75':False,'r90':False,
            'sb':False,'rb':False,'bxShape':False}
    for line in ff:
        errLine=line
        line=line.strip()
        if(line=='' or line[0]=='#'):
            continue
        if('=' not in line):
            print('Incorrect data format (equal sign missing) on line:\n'+errLine+'\nCalculations aborted')
            exit()
        vLine=line.split('=')
        vLine[0]=vLine[0].strip();vLine[1]=vLine[1].strip()
        if(vLine[0]=='' or vLine[1]=='' or len(vLine)!=2):
            print('incorrect(A) data format on line:\n'+errLine+'\nCalculations aborted')
            exit()        
        if(vLine[0] not in ['name','type','option']):            
            try:
                fVal=float(vLine[1])
                if(fVal<=0.0):
                    print('incorrect zero-value on line:\n'+errLine+'\nCalculations aborted')
                    exit()
                if(fVal>100.0):
                    print('unacceptable large value on line:\n'+errLine+'\nCalculations aborted')
                    exit()                
            except ValueError:
                if(vLine[1]=='*'):
                    fVal='*'
                else:    
                    print('incorrect(B) data format on line:\n'+errLine+'\nCalculations aborted')
                    exit()     
        if(vLine[0]=='name'):
            if(len(vLine[0])>36):
                print('Warning: name too long (only the first 30 chars are retained):\n')
            data['name']=vLine[1][0:30]; vCheck['name']=True;continue
        if(vLine[0]=='type'):
            if(vLine[1] in ['AL','FE']):
                data['type']=vLine[1]; vCheck['type']=True;continue
            else:
                print('incorrect material type on line:\n'+errLine+'\nCalculations aborted')
                exit()                    
        if(vLine[0]=='s0'):
            data['biax']['s0']=fVal; vCheck['s0']=True;continue
        if(vLine[0]=='r0'):
            data['biax']['r0']=fVal; vCheck['r0']=True;continue
        if(vLine[0]=='s90'):
            data['biax']['s90']=fVal; vCheck['s90']=True;continue 
        if(vLine[0]=='r90'):
            data['biax']['r90']=fVal; vCheck['r90']=True;continue
        if(vLine[0]=='sb'):
            data['biax']['sb']=fVal; vCheck['sb']=True;continue
        if(vLine[0]=='rb'):
            data['biax']['rb']=fVal; vCheck['rb']=True;continue
        if(vLine[0]=='bxShape'):
            data['biax']['bxShape']=fVal; vCheck['bxShape']=True;continue            
        if(vLine[0]=='s15'):
            data['UTstress'][0]=fVal; vCheck['s15']=True;continue   
        if(vLine[0]=='s30'):
            data['UTstress'][1]=fVal; vCheck['s30']=True;continue
        if(vLine[0]=='s45'):
            data['UTstress'][2]=fVal; vCheck['s45']=True;continue
        if(vLine[0]=='s60'):
            data['UTstress'][3]=fVal; vCheck['s60']=True;continue
        if(vLine[0]=='s75'):
            data['UTstress'][4]=fVal; vCheck['s75']=True;continue
        if(vLine[0]=='r15'):
            data['UTrvalue'][0]=fVal; vCheck['r15']=True;continue   
        if(vLine[0]=='r30'):
            data['UTrvalue'][1]=fVal; vCheck['r30']=True;continue
        if(vLine[0]=='r45'):
            data['UTrvalue'][2]=fVal; vCheck['r45']=True;continue
        if(vLine[0]=='r60'):
            data['UTrvalue'][3]=fVal; vCheck['r60']=True;continue
        if(vLine[0]=='r75'):
            data['UTrvalue'][4]=fVal; vCheck['r75']=True;continue 
        if(vLine[0]=='option'):
            if(vLine[1] not in ['0','1','2']):
                print('incorrect option (must be 0, 1 or 2) on line:\n'+errLine+'\nCalculations aborted')
                exit()
            data['option']=vLine[1]; vCheck['option']=True                
    for item in vCheck:
        if(not vCheck[item]):
            print(item+': no data provided\nCalculations aborted')
            exit()                       
                            
    

if __name__ == "__main__": 
    #data=ISO 
    #data=AA2090T3
    #data=DP980
    #data=AA3104H19
    #data=AISI304
    readData('mat000File.txt')
    print('----------------------------{}'.format(data['name']))
    [missingUAXdata,B90min,T45,B90max]=dataCheck()  
    ####bxParamPlot()  
    vBX=bxParam()
    Poly2()
    print('Smin={}, Smax={}, A4min={}, A4max={}'.format(vBX[0],vBX[1],vBX[2],vBX[3]))
    addedData=[]; strOption=''
    if(not missingUAXdata):
        vUAX=uaxParam()
    else:
        option=data['option']
        strOption='_opt_'+data['option']
        if(option=='0'):
            vUAX=uaxParamMissing()  ##with no added data
        if(option=='1'):            
            ##addedData format: list of (angle,yStress,rValue)
            midStressA=0.5*(data['biax']['s0']+data['UTstress'][2])
            midStressB=0.5*(data['UTstress'][2]+data['biax']['s90'])
            midRvalA=0.5*(data['biax']['r0']+data['UTrvalue'][2])
            midRvalB=0.5*(data['UTrvalue'][2]+data['biax']['r90'])
            addedData=[(PI4/2,midStressA,midRvalA),(3*PI4/2,midStressB,midRvalB)]
            vUAX=uaxParamAddedData(addedData)
        if(option=='2'):
            ##addedData format: list of (angle,yStress,rValue)
            midStressA=0.5*(data['biax']['s0']+data['UTstress'][2])
            midStressB=0.5*(data['UTstress'][2]+data['biax']['s90'])
            midRvalA=0.5*(data['biax']['r0']+data['UTrvalue'][2])
            midRvalB=0.5*(data['UTrvalue'][2]+data['biax']['r90'])
            addedData=[(PI4/2,midStressA,midRvalA),(3*PI4/2,midStressB,midRvalB)]
            LBD=0.3            
            s15=(1-LBD)*data['biax']['s0']+LBD*midStressA
            s30=(1-LBD)*midStressA+LBD*data['UTstress'][2]
            s60=(1-LBD)*midStressB+LBD*data['UTstress'][2]
            s75=(1-LBD)*data['biax']['s90']+LBD*midStressB
            r15=(1-LBD)*data['biax']['r0']+LBD*midRvalA
            r30=(1-LBD)*midRvalA+LBD*data['UTrvalue'][2]
            r60=(1-LBD)*midRvalB+LBD*data['UTrvalue'][2]
            r75=(1-LBD)*data['biax']['r90']+LBD*midRvalB            
            addedData=[(PI12,s15,r15),(2*PI12,s30,r30),(4*PI12,s60,r60),(5*PI12,s75,r75)]           
            vUAX=uaxParamAddedData(addedData)           
    print('convex={}, LL2={}, LL1={}, Nbisect={}'.format(vUAX[0],vUAX[1],vUAX[2],vUAX[3]))
    plotPoly6(missingUAXdata,addedData,strOption)
    with open('{name}_data{opt}.txt'.format(name=data['name'],opt=strOption),'w') as ff:
        ff.write('Hill48: min90 = {vv1:.4f}\nHill48: T45 = {vv2:.4f}\nHill48: max90 = {vv3:.4f}\n'.format(vv1=B90min,vv2=T45,vv3=B90max))
        ff.write('Smin={}, Smax={}, A4min={}, A4max={}\n'.format(vBX[0],vBX[1],vBX[2],vBX[3]))    
        ff.write('convex={}, LL2={}, LL1={}, Nbisect={}\n\n'.format(vUAX[0],vUAX[1],vUAX[2],vUAX[3]))
        for kk in range(1,17):
            print('a[{idx}] = {val:.4f}'.format(idx=kk,val=vCoeff[kk]))
            ff.write('a[{idx}] = {val:.4f}\n'.format(idx=kk,val=vCoeff[kk]))
            

    


