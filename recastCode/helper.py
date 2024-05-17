#!/usr/bin/env python3

import numpy as np
from scipy.interpolate import CubicSpline
import os


data_Iter5 = np.genfromtxt(os.path.join(os.path.dirname(__file__),'Iter_5.csv'),delimiter=',')
trackEff_F = CubicSpline(data_Iter5[:,0],data_Iter5[:,1],extrapolate=True)



def trackEff(track):

    r_cm = np.sqrt(track.X**2 + track.Y**2)/10.0
    # Extrapolate with a fixed eff
    if r_cm > 60.0:
        return 0.15
    eff = max(0.,trackEff_F(r_cm))
    eff = min(eff,1.0)
    return eff


def getJetTracks(jet,allTracks):

    jet_tracks = []
    for track in allTracks:
        if track.PT < 1.0:
            continue
        # Eta and Phi at the outer edge of the tracker:
        eta = track.EtaOuter
        phi = track.PhiOuter        
        dR = np.sqrt((jet.Eta-eta)**2 + (jet.Phi-phi)**2)
        if dR > 0.4:
            continue
        eff = trackEff(track)
        if np.random.uniform() > eff:
            continue
        ## Remove tracks with no hits in the pixel layers:
        # x = track.X
        # y = track.Y
        # vTrack = np.array([x,y])
        # rprod = np.linalg.norm(vTrack)
        # if rprod > 102.0:
            # continue
        jet_tracks.append(track)

    return jet_tracks

def getSigmaD0(track):

    x = track.X
    y = track.Y
    phi = track.Phi
    pT = track.PT
    vTrack = np.array([x,y])
    pTrack = np.array([pT*np.cos(phi),pT*np.sin(phi)])

    ## Assuming nominal values quoted in 1405.6569v2 (pg.21)
    a = 30e-3
    b = 100e-3
    d0Err = np.sqrt(a**2 + (b/track.PT)**2)
    # return d0Err

    ## Assuming the d0 error comes from the phi, and x,y resolutions:
    d0Err = np.sqrt(getSigmaXYZ(track)**2 + ((np.dot(vTrack,pTrack)/pT)**2)*getSigmaPhi(track)**2)
    return d0Err


def getSigmaXYZ(track):

    # Apply a resolution depending on the closest Pixel Layer
    # (could only find resolutions for Layers 1 and 3)
    # Note: 1405.6569v2 quotes a resolution of 10e-3 mm in the transverse coordinate and 20e-3-40e-3 mm in the longitudinal one
    vTrack = np.array([track.X,track.Y])
    rprod = np.linalg.norm(vTrack)
    if rprod < 44.0:
        return 33.92*1e-3 # Resolution for Barrel Pixel Layer 1 (according to 1710.03842 Fig.7)
    elif rprod < 102.0:
        return 12.49*1e-3 # Resolution for Barrel Pixel Layer 3 (according to 1710.03842 Fig.8)
    elif rprod < 550.0:
        return 38.0*1e-3 # Resolution for Tracker Inner Barrel (TIB)(according to 1405.6569v2 pg.4)
    else:
        return 47.0*1e-3 # Resolution for Tracker Outer Barrel (TOB)(according to 1405.6569v2 pg.4)    

def getSigmaPhi(track):
    
    return 15e-3 # 1e-4--1e-2 rad from 1405.6569 (Fig.14)

def getIP2D(tracks,smear=True):

    ipList = []    
    for track in tracks:
        d0Err = getSigmaD0(track)      

        ## Method 1: directly smear d0 by its estimated error
        # d0 = track.D0
        # if smear:
        #     d0 = np.random.normal(loc=track.D0,scale=d0Err)
        # if d0 == 0.0: # Hack to deal with zero impact parameter
        #     ipT = 0.0
        # else:
        #     ipT = np.log10(abs(d0)/d0Err)
        # ipList.append(ipT)
        
        ## Method 2: smear d0 by smearing the the track vertex and phi
        x = track.X
        y = track.Y
        phi = track.Phi
        pT = track.PT
        if smear: # PT has already been smeared
            x = np.random.normal(loc=x,scale=getSigmaXYZ(track))
            y = np.random.normal(loc=y,scale=getSigmaXYZ(track))
            phi = np.random.normal(loc=phi,scale=getSigmaPhi(track))

        vTrack = np.array([x,y])
        pTrack = np.array([pT*np.cos(phi),pT*np.sin(phi)])
        d0 = np.linalg.norm(np.cross(vTrack,pTrack))/pT        
        if d0 == 0.0: # Hack to deal with zero impact parameter
            ipT = 0.0
        else:
            ipT = np.log10(abs(d0)/d0Err)
        ipList.append(ipT)
    
    return np.median(ipList)

def getTheta2D(tracks,smear=True):

    thetaList = [] 
    # radial locations of first silicon pixels at CMS
    pxd = np.array([44,73,102])
    pxd.sort()
    for track in tracks:
        x = track.X
        y = track.Y
        phi = track.Phi
        pT = track.PT
        
        vTrack = np.array([x,y])
        rprod = np.linalg.norm(vTrack)
        pTrack = np.array([pT*np.cos(phi),pT*np.sin(phi)])
        # Check if the track starts after the first pixels or not
        if rprod > max(pxd):
            Rmin = rprod
        else:
            Rmin = pxd[np.argmax(pxd > rprod)]
        if rprod < Rmin: # production vertex has to be larger than inner most layer
            if rprod > 0.0:
                ctheta = np.dot(vTrack,pTrack)/(rprod*pT)
                # Compute track length required to hit first inner layer:
                l = np.sqrt(Rmin**2 - (1-ctheta**2)*rprod**2) - rprod*ctheta
            else:
                l = Rmin
            # Set the track production position to the inner most layer
            x = x + l*np.cos(track.Phi)
            y = y + l*np.sin(track.Phi)
        
        if smear: # PT has already been smeared
            x = np.random.normal(loc=x,scale=getSigmaXYZ(track))
            y = np.random.normal(loc=y,scale=getSigmaXYZ(track))
            phi = np.random.normal(loc=phi,scale=getSigmaPhi(track))

        # Now re-compute the angle:
        vTrack = np.array([x,y])
        pTrack = np.array([pT*np.cos(phi),pT*np.sin(phi)])
        rprod = np.linalg.norm(vTrack)        
        ctheta = np.dot(vTrack,pTrack)/(rprod*pT)
        ctheta = np.round(ctheta,6) # To avoid numerical instabilities
        theta = np.arccos(ctheta)
        ltheta = np.log10(theta)
        thetaList.append(ltheta)
        
    return np.median(thetaList)

def getAlpha(tracks,smear=True):

    sumPV = 0.0
    sumAll = 0.0
    for track in tracks:
        # if np.linalg.norm([track.X,track.Y,track.Z]) < 1e-5:
        sumAll += track.PT
        x = track.X
        y = track.Y
        z = track.Z
        if smear: # PT has already been smeared
            x = np.random.normal(loc=x,scale=getSigmaXYZ(track))
            y = np.random.normal(loc=y,scale=getSigmaXYZ(track))
            z = np.random.normal(loc=z,scale=getSigmaXYZ(track))

        vTrack = np.array([x,y,z])
        rprod = np.linalg.norm(vTrack)
        # if track.D0 < 0.0001 and rprod < 0.001:
        #     sumPV += track.PT
        if rprod < 1.0:
            sumPV += track.PT

    return sumPV/sumAll



