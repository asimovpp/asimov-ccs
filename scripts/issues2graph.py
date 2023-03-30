#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import csv
from pyvis.network import Network


HTML_EXPORT = "issue_graph.html"

hex_palette = ['#fab0e4', '#ffb482', '#8de5a1', '#ff9f9b', '#d0bbff', '#debb9b', '#cfcfcf', '#fffea3', '#b9f2f0', '#a1c9f4']
label_color = {"todo": "#E19F40",
               "perpetual":"#8600C4",
               "doing": "#4FAA4F" }

if len(sys.argv) > 1:
    csv_path = sys.argv[1]
else:
    csv_path = "asimov.csv"

issues = []

print("Reading {}".format(csv_path))

with open(csv_path) as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count = 0
    for row in csv_reader:
        if line_count != 0:
            issues.append({"ID"                : int(row[0]) ,  
                "URL"               : row[1] ,
                "Title"             : row[2] ,
                "State"             : row[3] ,
                "Description"       : row[4] ,
                "Author"            : row[5] ,
                "Author_Username"   : row[6] ,
                "Assignee"          : row[7] ,
                "Assignee_Username" : row[8] ,
                "Confidential"      : row[9] ,
                "Locked"            : row[10] , 
                "Due_Date"          : row[11] , 
                "Created_At"        : row[12] , 
                "Updated_At"        : row[13] , 
                "Closed_At"         : row[14] , 
                "Milestone"         : row[15] , 
                "Weight"            : row[16] , 
                "Labels"            : row[17].split(",") , 
                "Time_Estimate"     : row[18] , 
                "Time_Spent"        : row[19] })
        line_count += 1


assignee_list = []
for issue in issues:
    assignee = issue["Assignee"]
    if not assignee in assignee_list:
        assignee_list.append(assignee)



g = Network(width="100%", directed=True)

for issue in issues:
    color = "#2980B9"

    for label in issue["Labels"]:
        if label in label_color:
            color = label_color[label]
            break
        
    assignee_id = assignee_list.index(issue["Assignee"])
    g.add_node(issue["ID"], 
               title = issue["Assignee"],
               label = issue["Title"],
               color = color)

for issue in issues:

    if "DEPENDENCIES:" in issue["Description"]:
        for line in issue["Description"].split('\n'):
            dependencies = []
            if 'DEPENDENCIES:' in line:
                dependencies = line.split(':')[1].split()
                dependencies = [int(dep[1:]) for dep in dependencies] # remove the leading "#"
                break
            
        for dependency in dependencies:
            g.add_edge(issue["ID"], dependency, value = 3)
            print("{} -> {}".format(issue["ID"], dependency))


    # TODO: detect if edge already exists before adding it
    #if "REQUIRES:" in issue["Description"]:
    #    for line in issue["Description"].split('\n'):
    #        requires = []
    #        if 'REQUIRES:' in line:
    #            requires = line.split(':')[1].split()
    #            requires = [int(dep[1:]) for dep in dependencies] # remove the leading "#"
    #            break
    #        
    #    for required in requires:
    #        g.add_edge(required, issue["ID"])
    #        print("{} -> {}".format(required, issue["ID"]))




g.show(HTML_EXPORT)
print("Exporting to {}".format(HTML_EXPORT))


