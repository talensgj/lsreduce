# -*- coding: utf-8 -*-
"""
Created on Fri Dec  1 17:58:17 2017

@author: talens
"""

import time
import psutil
import logging
import traceback

WD = r'C:\Users\mascara\Documents\lsreduce'
RED = r'scripts\lsred.py'

PYTHON = r'python'

MAIL_USER = 'mascara@strw.leidenuniv.nl'
MAIL_PASSWORD = 'XOpl@net'
MAIL_RECIPIENTS = ['talens@strw.leidenuniv.nl']

# psutil.ABOVE_NORMAL_PRIORITY_CLASS
# psutil.BELOW_NORMAL_PRIORITY_CLASS
# psutil.HIGH_PRIORITY_CLASS
# psutil.IDLE_PRIORITY_CLASS
# psutil.NORMAL_PRIORITY_CLASS
# psutil.REALTIME_PRIORITY_CLASS

NORMAL_PRIORITY_CLASS = psutil.NORMAL_PRIORITY_CLASS
ABOVE_NORMAL_PRIORITY_CLASS = psutil.ABOVE_NORMAL_PRIORITY_CLASS
BELOW_NORMAL_PRIORITY_CLASS = psutil.BELOW_NORMAL_PRIORITY_CLASS

VERYLOW_IOPRIORITY_CLASS = 0
LOW_IOPRIORITY_CLASS = 1
NORMAL_IOPRIORITY_CLASS = 2

def start_program(process):
    
    # Start the program.
    try:
        proc = psutil.Popen(process['cmdline'], cwd=process['wd'])
    except Exception as e:
        log.error('Caught {} while starting control'.format(e))
        return None

    # Set the priorities.
    try:
        proc.nice(process['priority'])
        proc.ionice(process['iopriority'])
    except Exception as e:
        log.error('Caught {} while setting process priority'.format(e))

    return proc


def send_email(program, info):
    
    import smtplib
    from email.MIMEMultipart import MIMEMultipart
    from email.MIMEText import MIMEText
    
    msg = MIMEMultipart()
    msg['Subject'] = 'Crash in La Silla - ' + program
    msg['From'] = MAIL_USER
    msg['To'] = ', '.join(MAIL_RECIPIENTS)          
    
    message = 'Hi All,\n\n{} has crashed\n{}\n\nPlease Help,\n La Silla'.format(program, info)
    text = MIMEText(message)
    msg.attach(text)

    try:
        
        smtpserver = smtplib.SMTP("smtp.strw.leidenuniv.nl", 587)
        smtpserver.ehlo()
        smtpserver.starttls()
        smtpserver.login(MAIL_USER, MAIL_PASSWORD)  
        smtpserver.sendmail(MAIL_USER, MAIL_RECIPIENTS, msg.as_string())
        smtpserver.close()
        
    except Exception as e:
        
        fstr = 'The following email was not sent because {}\n{}'
        log.error(fstr.format(e, message))

    return

def find_existing_process(name):
    """This is probably rickety"""
    
    for p in psutil.process_iter():
        
        cmdline = p.cmdline()
        
        for c in cmdline:
            if name in c:
                return p
            
    return None

def monitor_loop(processes, max_restarts=1, max_restart_time=3*3600):

    for p in processes:
        processes[p]['proc'] = None
        processes[p]['time'] = time.time()
        processes[p]['count'] = 0
        processes[p]['warning'] = 0

    for p in processes:
        
        processes[p]['proc'] = find_existing_process(p)
        
        if processes[p]['proc'] is None:
            processes[p]['proc'] = start_program(processes[p])
        else:
            print 'Recovering {}'.format(p)

    while True:

        for p in processes:

            proc = processes[p]['proc']
            
            if proc is None or not proc.is_running():
                
                if processes[p]['count'] < max_restarts:
                    
                    processes[p]['time'] = time.time()
                    processes[p]['count'] += 1
                    
                    send_email(p, 'restarting...')
                    processes[p]['proc'] = start_program(processes[p])
                    
                else:
                    
                    if not processes[p]['warning']:
                        send_email(p, 'restart limit reached')
                        processes[p]['warning'] = True

            if (time.time() - processes[p]['time']) > max_restart_time:
                processes[p]['count'] = 0 
                processes[p]['warning'] = False

        time.sleep(10)
        
    return

if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser(description='Initialize the reduction loop.')
    parser.add_argument('camera', type=str, choices=['LSN', 'LSE', 'LSS', 'LSW', 'LSC'],
                        help='the camera to initialize the reduction loop for')
    args = parser.parse_args()

    # Create logger.
    logging.basicConfig()
    log = logging.getLogger('Guardian')

    # Create process.
    processes = dict()
    processes[RED] = {'cmdline':[PYTHON, RED, args.camera], 'wd': WD, 'priority': BELOW_NORMAL_PRIORITY_CLASS, 'iopriority': LOW_IOPRIORITY_CLASS}

    try:
        monitor_loop(processes)
    except Exception as e:
        with open(r'd:\guardian.txt', 'w') as f:
            f.write('Guardian died: {}\n{}'.format(e, traceback.format_exc()))

