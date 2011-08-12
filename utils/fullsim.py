from simulations import corr21cm, lofar, foregroundsck, pointsource



class FullSim(Map3d):

    
    components = { 'signal' : corr21cm.Corr21m,
                   'synchrotron' : lofar.LofarGSDE,
                   #'pointsources' : pointsource.DiMatteo,
                   'galacticfreefree' : foregroundsck.GalacticFreeFree,
                   'extragalacticfreefree' : foreground.ExtraGalacticFreeFree
                   }

    
    def _create(self):

        
