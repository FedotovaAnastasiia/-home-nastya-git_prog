#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  1 13:26:34 2022

@author: nastya
"""
import tkinter as tk
import matplotlib.pylab as plt
import numpy as np
from matplotlib import rcParams
from matplotlib import patches
rcParams['font.size']='15'
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from astropy.io import fits
import matplotlib
matplotlib.use("TkAgg")
import os
import seaborn as sns
from matplotlib.colors import ListedColormap
cmap = ListedColormap(sns.color_palette())
import sunpy.cm
import sunpy
import sunpy.map
import sunpy.data.sample
import cv2
import matplotlib.ticker as ticker
import sqlite3
from time import perf_counter
import tkinter.font as fnt
import tkinter.filedialog as fd
import moviepy.editor as moviepy
from os import listdir



path = '/home/nastya/New_evens/catalog/2022-01-12/'
flux_temper_path = '/home/nastya/test_new2/' # папака в которую будут сохранены результаты программы
path_database = '/home/nastya/test_new2/DATABASE.db' # путь к директории для создания базы данных

class Makeentry(tk.Frame):
     
    def __init__(self,parent,w=None,indata='',font=None):
        if w==None:
            w=10
            
        tk.Frame.__init__(self,parent,height=20,width=w, bg="blue")
                
        self.pack_propagate(0)
        self.entry=tk.Entry(self,bg='#fff')
        if font!=None:
            self.entry.config(font=font)
        self.entry.pack(expand=1,fill=tk.BOTH)
        self.entry.bind('<Any-KeyRelease>',self.nokeypress)
        self.data=indata

    def setbg(self,bg):
        self.entry.config(bg=bg)

    def nokeypress(self,ev=None):
        self.entry.delete(0,tk.END)
        self.entry.insert(0,self.data)

    def refresh(self,dataentry):
        self.data=dataentry
        self.entry.delete(0,tk.END)
        self.entry.insert(0,self.data)

    def getvalue(self):
        return self.entry.get()

def myfunct(ev=None):
    print( appl.x_start, appl.y_start,appl.x_end, appl.y_end) 


class MyApp(tk.Tk):
    def __init__(self):
        tk.Tk.__init__(self)
        self.enttrace=Makeentry(self,font='15')
        
        self.but_panel=tk.Frame(self)
        
        self.btn_file = tk.Button(self.but_panel,  text="Выберите файл", command=self.choose_file, font = fnt.Font(family="Helvetica",
                size=15, weight="bold"),  width=20,  height=1,  highlightbackground = "lime", background='azure', fg='black')
        
        self.but_ok=tk.Button(self.but_panel,text='Готово', font = fnt.Font(family="Helvetica", size=15, weight="bold"),  width=20,  height=1,  highlightbackground = "lime",
                                      background='azure', fg='black', command=self.quit)
        
        self.but_cansel=tk.Button(self.but_panel,text='Закрыть', font = fnt.Font(family="Helvetica",
                size=15, weight="bold"),  width=20,  height=1,  highlightbackground = "lime",                         
                  background='azure', fg='black', command=self.destroy)
                     
        
        self.fig = plt.figure(figsize=(9,9))
        self.axes=self.fig.add_subplot(1,1,1)
        
        self.canv_fig=FigureCanvasTkAgg(self.fig,self)
        
        self.is_drawind_rectangle=0
        self.x_start=None
        self.y_start=None
        self.x_end=None
        self.y_end=None

        self.canv_fig._tkcanvas.grid(row=0,column=0)
        self.title('Выберите область')
        
        self.enttrace.grid(row=1,column=0,columnspan=2,sticky='w'+'e'+'n'+'s')
        self.but_panel.grid(row=0,column=1,padx=10, pady = 10)
        
        self.btn_file.pack(padx=60, pady=10)
        self.but_ok.pack(padx=60, pady=10)
        self.but_cansel.pack(padx=60, pady=10)
        self.columnconfigure(0,weight=1)
        self.rowconfigure(0,weight=1);
        self.rowconfigure(1,weight=0)
        self.mouse_move_event_id=self.canv_fig.mpl_connect('button_press_event',self.get_start_point)
        self.mouse_move_event_id=self.canv_fig.mpl_connect('button_release_event',self.draw_rectangle)
        self.mouse_move_event_id=self.canv_fig.mpl_connect('motion_notify_event',self.tracecmd)
        self.but_ok.bind("<Button-1>",myfunct)
        
        
    def choose_file(self):
        self.filename = fd.askopenfilename(title = "Select file",filetypes = (("FITS Files","*.fits"),))
        self.draw_data()    
        self.canv_fig.draw()
     
    def draw_data(self):
        plt.figure(1)
        plt.cla();
        plt.grid()
        self.f = fits.open(self.filename)
        self.header = self.f[0].header
        self.hdu = self.f[0]
        self.s = self.hdu.data
        plt.title(str(self.header['OBS-FREQ']) + ' '+ 'GHz' +' '+ self.header['DATE-OBS'][0:10] +' '+ self.header['DATE-OBS'][11:19])
        real_x= np.arange(-255,257,1) * self.header['CDELT1']
        real_y= np.arange(-255,257,1) * self.header['CDELT2']
        dx = (real_x[1]-real_x[0])/2.
        dy = (real_y[1]-real_y[0])/2.
        self.extent = [real_x[0]-dx, real_x[-1]+dx, real_y[0]-dy, real_y[-1]+dy]
        plt.imshow(self.s,cmap=plt.cm.hot,vmin=-24.84, vmax=28839, origin = 'lower', extent= self.extent)
        plt.grid(color = "w",linewidth = "0.8", linestyle = "-.")   
       
          
    def draw_rectangle(self,ev=None):
        self.x_end=ev.xdata
        self.y_end=ev.ydata
        self.draw_data()
        rect=patches.Rectangle([self.x_start, self.y_start], self.x_end-self.x_start, self.y_end-self.y_start, color = 'lime', fill=False, linewidth = "3.2")
              
        self.axes.add_artist(rect)
        self.canv_fig.draw()
        self.is_drawind_rectangle=0

    def get_start_point(self,ev=None):
        self.x_start=ev.xdata
        self.y_start=ev.ydata
        self.is_drawind_rectangle=1

    def tracecmd(self,ev=None):
        x=ev.xdata
        y=ev.ydata
        self.trace_point=[x,y]
        if x!=None and y!=None:
            trace='x: {0}; y: {1}'.format('%.3f'%x,'%.3f'%y)
            self.enttrace.refresh(trace)
        else:
            self.enttrace.refresh('')

        if self.is_drawind_rectangle:
            self.draw_data()
            rect=patches.Rectangle([self.x_start,self.y_start],x-self.x_start, y-self.y_start,color='lime',
                                   fill=False,linewidth="3.2")
            self.axes.add_artist(rect)
            self.canv_fig.draw()
   
    
def user_pro(x_start, y_start, x_end, y_end, n = 512, step = 30):
    

    Big_SIZE = 15
    plt.rc('font', size=Big_SIZE)
    plt.rc('axes', titlesize=25)
       
    pix = 5
    D = 1
    
    sky = np.zeros((n, n))
    
    def sort_fits(path):
        array = os.listdir(path) # получаем список папок из \NN
        sortetd_file_list = sorted(array)
        new_path = []
        path_all_new = []
            
        for folder in sortetd_file_list:
            path_all = path + folder + '/I'
            path_all_new.append(path_all)
            file_list = os.listdir(path_all)
            time_sorted_list = sorted(file_list)
            new_path.append(time_sorted_list)
            
        return new_path, path_all_new
    
    new_path, path_all_new = sort_fits(path)
           
    def create_name_folder(path_all_new, new_path):
        
        dirPath_new = []
        video_name_new = []
        image_folder_new = []
        frequency = []
        
        for i in range(len(path_all_new)):
            f = fits.open(path_all_new[i]+'/'+new_path[i][0])
            header = f[0].header
            dirPath = flux_temper_path + header['DATE-OBS'][0:10]+'_'+str(header['OBS-FREQ']) + '_GHz'
            dirPath_new.append(dirPath)
            image_folder = dirPath
            video_folder = flux_temper_path+'VIDEO_SRH/'
            video_name = (video_folder+'/SRH_{}.avi'.format(header['DATE-OBS'][0:10]+'_'+str(header['OBS-FREQ']) + '_GHz'))
            video_name_new.append(video_name)
            image_folder_new.append(image_folder)
            header = f[0].header
            hdu = f[0]
            frequency_new = hdu.header['OBS-FREQ']
            frequency.append(frequency_new)
            
        return dirPath_new, video_folder, header, video_name_new, image_folder_new, frequency       
           
    dirPath_new, video_folder_new, header, video_name_new, image_folder_new, frequency  = create_name_folder(path_all_new,new_path)
          
    
    def create_folder(dirPath_new):
        for j in range(len(dirPath_new)):
            if not os.path.isdir(dirPath_new[j]):
                print('Создается папка для картинок!')
                os.mkdir(dirPath_new[j])
        else:
            print('Папка уже создана!')
            return dirPath_new
        
        
    def create_folder_video(video_folder_new):
          if not os.path.isdir(video_folder_new):
              print('Создается папка для видео!')
              os.mkdir(video_folder_new)
          else:
              print('Папка уже создана!')
          return video_folder_new
          
        
    create_folder(dirPath_new)   
    create_folder_video(video_folder_new)
    
    
    def convertion_arcsec(x_start, x_end, y_start, y_end, header):
        x1 = int(n/2+x_start/header['CDELT1'])
        x2=int(n/2+x_end/header['CDELT1'])
        y1=int(y_start/header['CDELT2']+n/2)
        y2=int(y_end/header['CDELT2']+n/2)
        print(x1, x2, y1, y2)
        return  x1, x2, y1, y2
    
    x1, x2, y1, y2 = convertion_arcsec(x_start, x_end, y_start, y_end, header)
        
        
    def model_disk(n,pix):
        radius1 = int(940 / pix + 0.5)
        for k in range(n):
            for q in range(n):
                x = k - n / 2
                y = q - n / 2
                l = np.sqrt(x ** 2 + y ** 2)
                sky[k, q] = (np.pi / 2 - np.arctan(D * (l - radius1))) / np.pi
                if sky[k, q] < 0.1:
                    sky[k, q] = 1
                    if sky[k, q] > 0.9:
                        sky[k, q] = 1
    
        sky_new = np.zeros((n, n))
        for e in range(n):
            for e1 in range(n):
                if sky[e, e1] == 1:
                    sky_new[e, e1] = sky[e, e1]
        return sky_new
    
    
    def plot_images(new_path):
        smaxout=[]
        time = []
        date_new = []
        s_new = []
        fig = plt.figure(figsize = (10,10))  
        fig.subplots_adjust(left=0.13, bottom=0, right=1, top=1, wspace=1, hspace=-0.5)
       
   
        t_load1 = 0.0
        t_load2 = 0.0
        t_load3 = 0.0
        t_load4 = 0.0
        t_plot = 0.0
        t_save = 0.0
        
        calculated_model_disk = model_disk(n, pix)[y1:y2, x1:x2]
  
        for i in range(len(path_all_new)):
            smax_out = []
            time_new = []
            #for j in range(len(new_path[i])):
            for j in range(len(new_path[i][0:5])):
                t_start = perf_counter()
                f = fits.open(path_all_new[i]+'/'+new_path[i][j])
                img_path = path_all_new[i]+'/'+new_path[i][j] 
                header = f[0].header
                hdu = f[0]
                date = (f[0].header['DATE-OBS'][11:19])
                date_all = f[0].header['DATE-OBS'][0:10]
                date_new.append(date_all)
                time_new.append(date)
                t_load1 += perf_counter() - t_start
                t_start = perf_counter()
                s = hdu.data[y1:y2, x1:x2]
                t_load2 += perf_counter() - t_start
                t_start = perf_counter()
                c1_old = s * calculated_model_disk
                t_load3 += perf_counter() - t_start
                t_start = perf_counter()
                s_new.append(c1_old)
                smax_out.append(c1_old.max())
                f.close()
                t_load4 += perf_counter() - t_start
                t_start = perf_counter()
                ax1 = fig.add_subplot(1,1,1)                                           
                smap = sunpy.map.Map(img_path)
                rect = patches.Rectangle([x_start, y_start], x_end-x_start, y_end-y_start, color = 'lime', fill=False, linewidth = "3.2")
                ax1.add_artist(rect)
                smap.plot_settings['title']= (str(header['OBS-FREQ']) + ' '+ 'GHz' +' '+ header['DATE-OBS'][0:10] +' '+ header['DATE-OBS'][11:19])
                smap.plot(cmap=plt.cm.hot, vmin=-24.84, vmax=28839)
                t_plot += perf_counter() - t_start
                t_start = perf_counter()
                plt.colorbar(shrink=0.71)
                plt.grid(color = "w",linewidth = "0.8", linestyle = "-.")
                             
                fig.savefig(dirPath_new[i]+'/SRH_{}.png'.format(header['DATE-OBS']+' '+ str(header['OBS-FREQ'])+' '+ 'GHz'), dpi =60.2, bbox_inches = 'tight',
                pad_inches = 0.1)
                
                
                plt.clf()
                t_save += perf_counter() - t_start
                
                          
                print(f'Opens: {t_load1:.1f} s, hdu.data[y1:y2, x1:x2]: {t_load2:.1f} s,'
                      f' model_disk: {t_load3:.1f} s, close(): {t_load4:.1f} s, plot: {t_plot:.1f} s, save: {t_save:.1f} s')
               
            smaxout.append(smax_out)
            time.append(time_new)
    
                                    
        return smaxout, time, s_new, date_new
             
    smaxout,time, s_new, date_new = plot_images(new_path)
     
     
    def Brightness_temperature(frequency):
        
        for freq in range(len(frequency)):
                         
            fig, (ax1) = plt.subplots(1, figsize=(12, 3))
            fig.subplots_adjust(hspace=0.6)
            ax1.plot(time[freq], smaxout[freq], 'r', label=str(frequency[freq]) + ' GHz')
            ax1.grid()
            ax1.xaxis.set_major_locator(ticker.MultipleLocator(step))
            ax1.set_title('SRH', fontweight="bold")
            ax1.set_ylabel('$Brightness$ $temperature$ $[K]$')
            ax1.set_xlabel('$Time$ ' ' $[UT]$')
            legend = ax1.legend()
            name=str.replace((str(header['DATE-OBS'])+' '+str(frequency[freq])+' '+'GHz'),':','_')
            fig.savefig(flux_temper_path+'/SRH_{0}.png'.format(name))
            plt.cla()
            plt.close()
                
    Brightness_temperature(frequency)
        
    
    def sort_fits_new(path):
        file_list = os.listdir(path)
        full_list = [os.path.join(path, i) for i in file_list]
        time_sorted_list = sorted(full_list)
        return time_sorted_list
    
    new_path = sort_fits_new(path)
    
    
    def create_video(image_folder_new):
        for i in range(len(image_folder_new)):
            images = [img for img in sort_fits_new(image_folder_new[i]) if img.endswith(".png")]
            frame = cv2.imread(os.path.join(image_folder_new[i], images[0]))
            height, width, layers = frame.shape
            video = cv2.VideoWriter(video_name_new[i], 0, 5, (width,height)) # 5 - скорость видео
            for image in images:
                video.write(cv2.imread(os.path.join(image_folder_new[i], image)))  
                cv2.destroyAllWindows()
            video.release()    
    create_video(image_folder_new)
            
    
    def zip_video(video_folder_new):
        file_list = os.listdir(video_folder_new)
        for i in range(len(file_list)):
            clip = moviepy.VideoFileClip(video_folder_new + file_list[i])
            clip.write_videofile(video_folder_new+(str(file_list[i][:-4])+'.mp4'))
              
    
    zip_video(video_folder_new)
    
    def delite_avi(video_folder_new):
        for file_name in listdir(video_folder_new):
            if file_name.endswith('.avi'):
                os.remove(video_folder_new + file_name)
     
    delite_avi(video_folder_new)
    
    
    def General_Base(path_database):   
        
        connection = sqlite3.connect(path_database, timeout=1)
        cursor = connection.cursor()
        tb_create = ('''CREATE TABLE events
                          (date int, time_start int, time_end int, arc_center_x int, arc_center_y int, width int, height int )''')
                         
        tb_exists = ("SELECT name FROM sqlite_master WHERE type='table' AND name='events'")
        
        if not connection.execute(tb_exists).fetchone():
            connection.execute(tb_create)
            
        date = date_new[0]
        time_start = time[0][0]
        time_end = time[-1][-1]
        arc_center_x = (x_start+x_end)/2
        arc_center_y = (y_start+y_end)/2
        width = abs(x_start-x_end)
        height = abs(y_start-y_end)
        
        cursor.execute("INSERT INTO events VALUES ( ?, ?, ?, ?, ?, ?, ?)", (date, time_start, time_end, arc_center_x, arc_center_y, width, height))
        cursor.execute("SELECT * FROM events")
        rows = cursor.fetchall()
        
        for row in rows:
         	print(row)
        
        tb_create = ('''CREATE TABLE observations
                          (frequency real, br_temperatures int, path_video text, link_video text, link_images text)''')
                          
        tb_exists = ("SELECT name FROM sqlite_master WHERE type='table' AND name='observations'")
        
        if not connection.execute(tb_exists).fetchone():
            connection.execute(tb_create)
            
        for freq in range(len(frequency)):
            
            br_temperatures = float('{:.4f}'.format(max(smaxout[freq])))
            path_video = video_folder_new
            link_video = video_name_new[freq]
            link_images = dirPath_new[freq]
            cursor.execute("INSERT INTO observations VALUES (?, ?, ?, ?, ?)", (frequency[freq], br_temperatures, path_video, link_video, link_images))
            cursor.execute("SELECT * FROM observations ")
            rows = cursor.fetchall()
            
        connection.commit()
        cursor.close()
        connection.close()
        
        for row in rows:
            print(row)
    
    General_Base(path_database)     
    
    
if __name__ == "__main__":
     
    appl = MyApp()
    appl.mainloop()    
    
    user_pro(x_start=appl.x_start, y_start = appl.y_end, x_end = appl.x_end, y_end = appl.y_start, n = 512, step = 30)   
    
    

    