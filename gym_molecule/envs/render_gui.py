from tkinter import *
from PIL import ImageTk, Image


class Render():
    
    def __init__(self, master):
        self.mols = []
        self.obs = []
        self.frames = []
        
        self.index = -1
        self.current = -1

        self.topFrame = Frame(master)
        self.topFrame.pack()
        self.bottomFrame = Frame(master)
        self.bottomFrame.pack(side = BOTTOM)
        self.midFrame = Frame(master)
        self.midFrame.pack(side = BOTTOM)

        self.molLabel = Label(self.topFrame)
        self.molLabel.pack()
        self.graph = Label(self.midFrame)
        self.graph.pack()

    def update(self, mol):
        self.mols.append(mol)
        self.index += 1

        for i in range(360):
            pic = "/gym_molecule/envs/resources/"+str(i)+".png"
            image = Image.open(pic).resize((400,400), Image.ANTIALIAS)
            frame = ImageTk.PhotoImage(image)
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

            self.frames = self.obs[self.current]
            self.midFrame.after(0, self.updateGif, 0)

    def prev(self):
        if self.current > 0:
            self.current -= 1
            pic = ImageTk.PhotoImage(self.mols[self.current], master = self.molLabel)
            self.molLabel.configure(image=pic)
            self.molLabel.image = pic

            self.frames = self.obs[self.current]
            self.midFrame.after(0, self.updateGif, 0)

    def render(self):

        if len(self.mols) != 0:
            self.current = self.index
            pic = ImageTk.PhotoImage(self.mols[self.index])
            self.molLabel.configure(image=pic)
            self.molLabel.image = pic

            prevButton = Button(self.bottomFrame, text="PREVIOUS", command = self.prev)
            prevButton.pack(side=LEFT)
            nextButton = Button(self.bottomFrame, text="NEXT", command = self.nxt)
            nextButton.pack(side=RIGHT)
            self.midFrame.after(0, self.updateGif, 0)

    def reset(self):
        self.mols = []
        self.obs = []
        self.frames = []

