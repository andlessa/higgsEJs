#!/usr/bin/env python3

import numpy as np
from scipy.interpolate import CubicSpline
import os
import itertools
np.random.seed(0)


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


def getJetTracks(jet,allTracks,Rmax=None):

    jet_tracks = []
    for track in allTracks:
        if track.PT < 1.0:
            continue
        # Eta and Phi at the outer edge of the tracker:
        Lxy = np.sqrt(track.X**2 + track.Y**2)
        if Rmax is not None and Lxy > Rmax:
            continue
        eta = track.EtaOuter
        phi = track.PhiOuter        
        dR = np.sqrt((jet.Eta-eta)**2 + (jet.Phi-phi)**2)
        if dR > 0.4:
            continue
        eff = trackEff(track)
        if np.random.uniform() > eff:
            continue
        jet_tracks.append(track)

    return jet_tracks

def getSigmaD0(track,method):

    x = track.X
    y = track.Y
    phi = track.Phi
    pT = track.PT
    vTrack = np.array([x,y])
    pTrack = np.array([pT*np.cos(phi),pT*np.sin(phi)])

    ## Method A: Assuming nominal values quoted in 1405.6569v2 (pg.21)
    if method == 'smearD0':
        a = 30e-3
        b = 10e-3
        d0Err = np.sqrt(a**2 + (b/pT)**2)
    ##  Method B: Assuming the d0 error comes from the phi, and x,y resolutions:
    elif method == 'smearVtx':
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

def getD0(track,smear=True,method='smearVtx'):

    if not smear:
        d0 = abs(track.D0)
    elif method == 'smearD0':
    ## Method 1: directly smear d0 by its estimated error
        d0 = abs(track.D0)
        d0Err = getSigmaD0(track,method='smearD0')
        d0 = abs(np.random.normal(loc=track.D0,scale=d0Err))
    elif method == 'smearVtx':    
    ## Method 2: smear d0 by smearing the the track vertex and phi
        x = track.X
        y = track.Y
        phi = track.Phi
        pT = track.PT
        x = np.random.normal(loc=x,scale=getSigmaXYZ(track))
        y = np.random.normal(loc=y,scale=getSigmaXYZ(track))
        phi = np.random.normal(loc=phi,scale=getSigmaPhi(track))
        vTrack = np.array([x,y])
        pTrack = np.array([pT*np.cos(phi),pT*np.sin(phi)])
        d0 = np.linalg.norm(np.cross(vTrack,pTrack))/pT        

    return d0

def getIP2D(tracks,smear=True,method='smearVtx'):

    ipList = []    
    for track in tracks:
        d0Err = getSigmaD0(track,method=method)      
        d0 = getD0(track,smear=smear,method=method)
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

def getFourMom(pObj,mass=0.0):

    if hasattr(pObj,'Px'):
        return np.array([pObj.E,pObj.Px,pObj.Py,pObj.Pz])
    elif hasattr(pObj,'PT'):
        p3 = np.array([pObj.PT*np.cos(pObj.Phi),
                       pObj.PT*np.sin(pObj.Phi),
                       pObj.PT*np.sinh(pObj.Eta)])
        e = np.sqrt(mass**2 + np.dot(p3,p3))
        p4 = np.concatenate(([e],p3))
        return p4
    else:
        return None

def eventReader(tfileObj,nevts=-1,method='smearVtx',Rmax=None):
    """
    Reads a Delphes ROOT TFile object and collects the relevant distributions
    """

    tree = tfileObj.Get("Delphes")
    if nevts < 0:
        nevts = tree.GetEntries()

    
    vars = ['IP2D', 'theta2D', 'alpha','pT','d0','sigmaD0','weights']
    resDict = {v : [] for v in vars}
    resDict['method'] = method
    resDict['Rmax'] = Rmax

    for ievt in range(nevts):
        tree.GetEntry(ievt)   
        # weightPB = tree.Event.At(0).Weight/nevts
        weightPB = 1.0
        jets = tree.Jet
        # jets = tree.GenJet
        tracks = tree.Track        

        #Apply lepton cut requirement:
        electrons = tree.Electron
        muons = tree.Muon
        passLepton = False
        if len(electrons) >= 2:
            ePlus = [e for e in electrons if e.Charge > 0]
            eMinus = [e for e in electrons if e.Charge < 0]
            for eP,eM in itertools.product(ePlus,eMinus):
                p1 = getFourMom(eP,mass=0.511e-3)
                p2 = getFourMom(eM,mass=0.511e-3)
                p12 = p1+p2       
                mll = np.sqrt(p12[0]**2 - np.dot(p12[1:],p12[1:]))
                pTll = np.linalg.norm(p12[1:3])
                if pTll > 100.0 and (70.0 < mll < 110.):
                    passLepton = True
                    break
        if not passLepton and len(muons) >= 2:
            muPlus = [mu for mu in muons if mu.Charge > 0]
            muMinus = [mu for mu in muons if mu.Charge < 0]
            for muP,muM in itertools.product(muPlus,muMinus):
                p1 = getFourMom(muP,mass=106e-3)
                p2 = getFourMom(muM,mass=106e-3)
                p12 = p1+p2       
                mll = np.sqrt(p12[0]**2 - np.dot(p12[1:],p12[1:]))
                pTll = np.linalg.norm(p12[1:3])
                if pTll > 100.0 and (70.0 < mll < 110.):
                    passLepton = True
                    break
        
        if not passLepton:
            continue

        for j in jets:
            if j.PT < 35.0:
                continue
            if abs(j.Eta) > 2.4:
                continue
            jet_tracks = getJetTracks(j,tracks,Rmax=Rmax)
            if len(jet_tracks) == 0:
                continue
            for track in jet_tracks:
                resDict['d0'].append(getD0(track,smear=True,method=method))
                resDict['sigmaD0'].append(getSigmaD0(track,method=method))
            resDict['pT'].append(j.PT)
            resDict['alpha'].append(getAlpha(jet_tracks))
            resDict['IP2D'].append(getIP2D(jet_tracks,smear=True,method=method))
            resDict['theta2D'].append(getTheta2D(jet_tracks))
            resDict['weights'].append(weightPB)

    for v in vars:
        resDict[v] = np.array(resDict[v])

    return  resDict
    


