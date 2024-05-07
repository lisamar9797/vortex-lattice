"""
Date: Q4 2023 - Q1 2024
Author: Lisa Martinez
Institution: Technical University of Madrid

Description: This subroutine calculates the induced velocities Ux,Uy,Uz from a
semi infinitly vortex in the point px,py,pz. The calculations are made for 1
blade. The routine uses the helix radius r, pitch p and longitudinal starting
point x as input. The vortex starts at -infinity and stops at -x. The
calculations follows the procedure outlined by de Jong.
"""


import numpy as np


def De_Jong(x, r, p, phi, px, py, pz):

    from sources.Ci_P import Ci
    from sources.Si_P import si

    Ux = 0.0               # Initialization  of the variable U_x
    Uy = 0.0               # Initialization  of the variable U_y
    Uz = 0.0               # Initialization  of the variable U_z

    # CONSTANTS

    pbar=2*np.pi/p
    xbar=pbar*x
    x2bar=2*xbar
    xtld=pbar*px
    x2tld=2*xtld

    xsum=xbar+xtld
    xsum2=2*xsum*xsum
    xsum3=1.5*xsum2*xsum
    xsum4=xsum2*xsum2

    x2sum=2*xsum
    x2sum2=2*x2sum*x2sum
    x2sum3=1.5*x2sum2*x2sum
    x2sum4=x2sum2*x2sum2

    r2=r*r
    py2=py*py
    pz2=pz*pz
    rbar=r2+py2+pz2
    p2=pbar*pbar
    p3r2=p2*r2*pbar
    p4r2=3*p3r2*pbar
    p5r2=p4r2*pbar
    p2x2=p2/xsum2
    p4rb=3*p2*p2*rbar/xsum4

    cxtld=np.cos(xtld)
    sxtld=np.sin(xtld)
    cx2tld=np.cos(x2tld)
    sx2tld=np.sin(x2tld)
    cxbar=np.cos(xbar)
    sxbar=np.sin(xbar)
    cx2bar=np.cos(x2bar)
    sx2bar=np.sin(x2bar)
    cxbar2=cxbar*cxbar
    sxbar2=sxbar*sxbar

    cphi=np.cos(phi)
    sphi=np.sin(phi)
    cphi2=cphi*cphi
    sphi2=sphi*sphi
    phi2=2*phi
    c2phi=np.cos(phi2)
    s2phi=np.sin(phi2)

    Ci1 = Ci(xsum)
    Ci2 = Ci(x2sum)
    si1 = si(xsum)
    si2 = si(x2sum)

    cos11  = -cxtld*Ci1-sxtld*si1
    cos112 = -cx2tld*Ci2-sx2tld*si2

    sin11  = sxtld*Ci1-cxtld*si1
    sin112 = sx2tld*Ci2-cx2tld*si2

    cos12  = cxbar/xsum-sin11
    cos122 = cx2bar/x2sum-sin112

    sin12  = sxbar/xsum+cos11
    sin122 = sx2bar/x2sum+cos112

    cos13  = cxbar/xsum2-0.5*sin12
    cos132 = cx2bar/x2sum2-0.5*sin122

    sin13  = sxbar/xsum2+0.5*cos12
    sin132 = sx2bar/x2sum2+0.5*cos122

    cos23  = cxbar2/xsum2-sin122
    sin23  = sxbar2/xsum2+sin122

    cos14  = cxbar/xsum3-sin13/3
    cos142 = cx2bar/x2sum3-sin132/3

    sin14  = sxbar/xsum3+cos13/3
    sin142 = sx2bar/x2sum3+cos132/3

    cos24  = cxbar2/xsum3-4*sin132/3

    sin24  = sxbar2/xsum3+4*sin132/3

    cos15  = cxbar/xsum4-0.25*sin14
    cos152 = cx2bar/x2sum4-0.25*sin142

    sin15  = sxbar/xsum4+0.25*cos14
    sin152 = sx2bar/x2sum4+0.25*cos142

    cos25  = cxbar2/xsum4-2*sin142
    sin25  = sxbar2/xsum4+2*sin142

    # CALCULATION OF UX

    p3 = p2*pbar
    p4 = p3*pbar
    p5 = p4*pbar
    rbar3r = rbar+2*r2

    c0 = -p5r2*rbar/xsum4/2+p3r2/xsum2
    c1 = -p3*r*pz
    c2 = -p3*r*py
    c3 = 3*p5*r*pz*rbar3r/2
    c4 = 3*p5*r*py*rbar3r/2
    c5 = -p5r2*pz2
    c6 = -p5r2*py2
    c7 = -16*p5r2*py*pz

    Ux = (c0+(c1*cphi-c2*sphi)*cos13 + (c2*cphi+c1*sphi)*sin13 + (c3*cphi-c4*sphi)*cos15+(c4*cphi+c3*sphi)*sin15 +
          (c5*cphi2+c6*sphi2)*cos25 + (c6*cphi2+c5*sphi2)*sin25 + ((c5-c6)*8*s2phi+c7*c2phi)*sin152 - c7*s2phi*cos152)/4/np.pi

    # CALCULATION OF UY

    d0 =-p2*pz/xsum2+1.5*p4*pz*rbar/xsum4
    d1 = p2*r
    d2 = d1
    d3 = -3*p4*r*rbar/2
    d4 = 4*p4r2*pz
    d5 = p4r2*py
    d6 = -3*p4*r*py*pz
    d7 = 8*p4r2*py
    d8 = -3*p4*r*(py2+3*pz2+r2)/2
    d9 = p4r2*pz

    Uy = (d0+d1*cphi*cos13+d1*sphi*sin13+(d8*cphi-d6*sphi)*cos15 + (d8*sphi+d6*cphi)*sin15+0.25*d4*cphi*cphi*cos25 +
          0.25*d4*sphi*sphi*sin25 + (8*d5*c2phi+2*d4*s2phi)*sin152 - d7*s2phi*cos152-d1*sphi*cos12+d1*cphi*sin12-
          d3*sphi*cos14 + d3*cphi*sin14+d5*sphi*sphi*cos24+d5*cphi*cphi*sin24 + (d4*c2phi-0.5*d7*s2phi)*sin142-
          d4*s2phi*cos142)/4/np.pi

    # CALCULATION OF UZ

    e0 = p2*py/xsum2-1.5*p4*py*rbar/xsum4
    e1 = p2*r
    e2 = -e1
    e3 = -3*p4*r*rbar/2
    e4 = 4*p4r2*py
    e5 = p4r2*pz
    e6 = 3*p4*r*py*pz
    e7 = 3*p4*r*(3*py2+pz2+r2)/2
    e8 = -8*p4r2*pz
    e9 = -p4r2*py

    Uz = (e0+e1*cphi*cos12+e1*sphi*sin12+e1*sphi*cos13-e1*cphi*sin13 + e3*cphi*cos14+e3*sphi*sin14+
         (e4*c2phi-0.5*e8*s2phi)*sin142 - e4*s2phi*cos142+e5*cphi*cphi*cos24+e5*sphi*sphi*sin24 -
          e7*sphi*cos15+e7*cphi*sin15+e6*cphi*cos15+e6*sphi*sin15 - e8*s2phi*cos152+
          (2*e4*s2phi+e8*c2phi)*sin152 + e9*sphi*sphi*cos25+e9*cphi*cphi*sin25)/4/np.pi

    return (Ux,Uy,Uz)
