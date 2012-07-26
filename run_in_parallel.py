#!/usr/local/bin/python

'''
#
# run_in_parallel.py
#
# Purpose:
#       to run one or more reports in parallel, where those reports are
#	specified on the command-line
#
# Usage:
#       run_in_parallel.py <script 1> [<script 2>...]
#
# Example:
#	run_in_parallel.py weekly/*.py
#
# History:
# 	jsb : 7-23-12 : initial creation
'''

import sys
import os
import glob
import time
import getopt
import Dispatcher

###--- Globals ---###

USAGE = '''Usage: %s [-p <process count>] <script 1> [<script 2>...]
    Purpose:
	Runs one or more report scripts in parallel to reduce overall runtime.
    Arguments:
	-p : number of parallel processes to use (default = 2)
''' % sys.argv[0]

REPORTS = []				# names of reports to be run
START_TIME = time.time()		# time (in ms) at which we began
IDS = []				# list of (report ID, report name)
DISPATCHER = Dispatcher.Dispatcher(2)	# runs the reports in parallel

###--- Functions ---###

def bailout (msg, showUsage = True):
	if showUsage:
		print USAGE,
	print 'Error: %s' % msg
	sys.exit(1)

def elapsedTime():
	return '%8.3f' % (time.time() - START_TIME)

def reportProgress():
	global IDS, DISPATCHER

	remaining = []		# reports not yet completed: (id, report name)

	for (report_id, report) in IDS:
		if DISPATCHER.getStatus(report_id) != Dispatcher.FINISHED:
			remaining.append ( (report_id, report) )
		else:
			# this one is finished; if it exited abnormally, then
			# we need to bail out

			if DISPATCHER.getReturnCode(report_id):
				DISPATCHER.terminateProcesses()
				bailout ('Error: %s failed; terminating other reports' % report)
			print '%s sec : Finished %s in %0.3f' % (
				elapsedTime(), report,
				DISPATCHER.getElapsedTime(report_id) )

	IDS = remaining
	return

###--- Main ---###

# process command-line

try:
	(options, args) = getopt.getopt (sys.argv[1:], 'p:')
except getopt.GetoptError:
	bailout ('Improper command-line')

if len(args) < 1:
	bailout ('Too few arguments')

for (option, value) in options:
	if option == '-p':
		DISPATCHER = Dispatcher.Dispatcher (int(value))
	else:
		bailout ('Unknown flag: %s' % option)

# assemble our list of reports from the command-line arguments

for report in args:
	# if this is the name of a report, then just add it
	if os.path.isfile(report):
		REPORTS.append(report)

	# if this is the path to a directory, then add all of its contained
	# reports (for convenience)
	elif os.path.isdir(report):
		for script in glob.glob (os.path.join (report, '*.py')):
			REPORTS.append (script)

	else:
		# this may have a wildcard contained in it, which we need to
		# evaluate and add individual matching items
		for script in glob.glob (report):
			REPORTS.append (script)

# schedule the reports to run in parallel

IDS = []

cmd = './runone_report.csh %s'

for report in REPORTS:
	report_id = DISPATCHER.schedule(cmd % report)
	IDS.append ( (report_id, report) )

DISPATCHER.wait(reportProgress)

