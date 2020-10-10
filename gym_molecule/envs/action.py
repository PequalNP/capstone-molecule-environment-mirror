class Action():
    def __init__(self):
        self.action_c = ''
        self.pos = ''   #front or back
        self.mol = ''
        self.query = ''
        self.isSmarts  = False

    def setAction(self,action,pos='front',query='',mol='C',isSmarts=False): #mol
        self.action_c  = action
        self.mol       = mol
        self.pos       = pos 
        self.query     = query
        self.isSmarts  = isSmarts       