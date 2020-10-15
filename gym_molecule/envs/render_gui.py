from tkinter import *
from PIL import ImageTk, Image


class Render():
    
    def __init__(self, master):
        self.mols = []
        self.obs = []
        self.frames = []
        
        self.index = -1
        self.current = -1

#         self.root = Tk()
        self.topFrame = Frame(master)
        self.topFrame.pack()
        self.bottomFrame = Frame(master)
        self.bottomFrame.pack(side = BOTTOM)
        self.midFrame = Frame(master)
        self.midFrame.pack(side = BOTTOM)

        self.molLabel = Label(self.topFrame)
        self.molLabel.pack()
        # self.graphL = Label(self.midFrame)
        # self.graphL.pack(side=LEFT, fill=X)
        self.graph = Label(self.midFrame)
        self.graph.pack()

    def update(self, mol):
        self.mols.append(mol)
        # self.obs1.append(graph1)
        # self.obs2.append(graph2)
        self.index += 1
        # self.frames = [PhotoImage(file='./resources/regr.gif',format = 'gif -index %i' %(i)) for i in range(360)]
        for i in range(360):
            pic = "./resources/"+str(i)+".png"
            frame = ImageTk.PhotoImage(Image.open(pic))
            self.frames.append(frame)
        self.obs.append(self.frames)

    def updateGif(self, ind):
        frame = self.frames[ind]
        if ind == 359:
            ind = 0
        else:
            ind +=1
        self.graph.configure(image = frame)
        self.graph.image = frame
        self.midFrame.after(40, self.updateGif, ind)

    def nxt(self):
        if self.current < self.index:
            self.current += 1
            pic = ImageTk.PhotoImage(self.mols[self.current])
            self.molLabel.configure(image=pic)
            self.molLabel.image = pic
            # pic = ImageTk.PhotoImage(self.obs1[self.current])
            # self.graphL.configure(image=pic)
            # self.graphL.image = pic
            # pic = ImageTk.PhotoImage(self.obs2[self.current])
            #  self.graph.configure(image=pic)
            # self.graph.image = pic
            self.frames = self.obs[self.current]
            self.midFrame.after(0, self.updateGif, 0)

    def prev(self):
        if self.current > 0:
            self.current -= 1
            pic = ImageTk.PhotoImage(self.mols[self.current], master = self.molLabel)
            self.molLabel.configure(image=pic)
            self.molLabel.image = pic
            # pic = ImageTk.PhotoImage(self.obs1[self.current])
            # self.graphL.configure(image=pic)
            # self.graphL.image = pic
            # pic = ImageTk.PhotoImage(self.obs2[self.current])
            # self.graph.configure(image=pic)
            # self.graph.image = pic
            self.frames = self.obs[self.current]
            self.midFrame.after(0, self.updateGif, 0)

    def render(self):
        # image = ImageTk.PhotoImage(img)
        # sprite = self.canvas.create_image(150,150, image= image)
        # img = Draw.MolToImage(self.current_molecule, size=(300,300), kekulize=True, wedgeBonds=True)
        # root = Tk()

        if len(self.mols) != 0:
            self.current = self.index
            pic = ImageTk.PhotoImage(self.mols[self.index])
            self.molLabel.configure(image=pic)
            self.molLabel.image = pic
            # pic = ImageTk.PhotoImage(self.obs1[self.index])
            # self.graphL.configure(image=pic)
            # self.graphL.image = pic
            # pic = ImageTk.PhotoImage(self.obs2[self.index])
            # self.graph.configure(image=pic)
            # self.graph.image = pic

            prevButton = Button(self.bottomFrame, text="PREVIOUS", command = self.prev)
            prevButton.pack(side=LEFT)
            nextButton = Button(self.bottomFrame, text="NEXT", command = self.nxt)
            nextButton.pack(side=RIGHT)
            self.midFrame.after(0, self.updateGif, 0)

    def reset(self):
        self.mols = []
        self.obs = []
        self.frames = []

