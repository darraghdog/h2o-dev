#!/usr/bin/python

import sys
import os
import requests
import urllib
import re


g_user = None
g_pass = None
g_sprint = None
g_csv = False
g_verbose = False


class Person:
    def __init__(self, name):
        self.name = name
        self.resolved_list = []
        self.resolved_story_points = 0
        self.unresolved_list = []
        self.unresolved_story_points = 0

    def add(self, issue):
        story_points = Person._get_story_points(issue)
        resolution = issue[u'fields'][u'resolution']
        if resolution is None:
            self.unresolved_story_points += story_points
            self.unresolved_list.append(issue)
        else:
            self.resolved_story_points += story_points
            self.resolved_list.append(issue)

    def emit_issues_csv(self):
        print (self.name + "," + str(len(self.resolved_list)) + "," + str(len(self.unresolved_list)))

    def emit_storypoints_csv(self):
        print (self.name + "," + str(self.resolved_story_points) + "," + str(self.unresolved_story_points))

    def emit_barchart(self):
        print("")
        print("-----" + self.name + "-----")
        Person._printbar("  resolved", self.resolved_story_points, "R")
        if (g_verbose):
            print("")
            for issue in self.resolved_list:
                story_points = Person._get_story_points(issue)
                print("{0:14s}{1:11s} ({2:.1f}): {3}".format("", issue[u'key'], story_points, issue[u'fields'][u'summary']))
        if (g_verbose):
            print("")
        Person._printbar("unresolved", self.unresolved_story_points, "U")
        if (g_verbose):
            print("")
            for issue in self.unresolved_list:
                story_points = Person._get_story_points(issue)
                print("{0:14s}{1:11s} ({2:.1f}): {3}".format("", issue[u'key'], story_points, issue[u'fields'][u'summary']))

    @staticmethod
    def _printbar(label, value, char):
        bars_per_day = 4
        sys.stdout.write(label + ":  ")
        num_bars = int(value * bars_per_day)
        i = 0
        while (num_bars > 0):
            if (i % bars_per_day == 0):
                sys.stdout.write(" ")
            sys.stdout.write(char)
            num_bars -= 1
            i += 1
        sys.stdout.write("\n")

    @staticmethod
    def _get_story_points(issue):
        story_points = issue[u'fields'][u'customfield_10004']
        if story_points is None:
            story_points = 0
        else:
            story_points = float(story_points)
        return story_points


class PeopleManager:
    def __init__(self):
        self.people_map = {}

    def add(self, issue):
        assignee = issue[u'fields'][u'assignee']
        if (assignee is None):
            print("ERROR: assignee is none for issue: " + str(issue))
            sys.exit(1)
        assignee_name = issue[u'fields'][u'assignee'][u'name']
        person = self.find(assignee_name)
        person.add(issue)

    def find(self, name):
        if name not in self.people_map:
            person = Person(name)
            self.people_map[name] = person
        person = self.people_map[name]
        return person

    def emit(self):
        global g_csv

        if g_csv:
            print("STORY_POINTS")
            print("name,resolved,unresolved")
            for key in sorted(self.people_map.keys()):
                person = self.people_map[key]
                person.emit_storypoints_csv()
            print("")
            print("NUMBER_OF_ISSUES")
            print("name,resolved,unresolved")
            for key in sorted(self.people_map.keys()):
                person = self.people_map[key]
                person.emit_issues_csv()
        else:
            for key in sorted(self.people_map.keys()):
                person = self.people_map[key]
                person.emit_barchart()


def usage():
    print("")
    print("usage:  " + g_script_name + " --user username --pass password --sprint sprintname")
    print("")
    sys.exit(1)


def unknown_arg(s):
    print("")
    print("ERROR: Unknown argument: " + s)
    print("")
    usage()


def parse_config_file():
    global g_user
    global g_pass
    global g_sprint
    home = os.path.expanduser("~")
    config_file = os.path.join(home, ".h2o_jira")
    if (os.path.exists(config_file)):
        with open(config_file) as f:
            for line in f:
                match_groups = re.search(r"\s*(\S+)\s*=\s*([\S\s]+)\s*", line.rstrip())
                if (match_groups is not None):
                    key = match_groups.group(1)
                    value = match_groups.group(2)
                    if (key == "user"):
                        g_user = value
                    elif (key == "pass"):
                        g_pass = value
                    elif (key == "sprint"):
                        g_sprint = value


def parse_args(argv):
    global g_user
    global g_pass
    global g_sprint
    global g_csv
    global g_verbose

    i = 1
    while (i < len(argv)):
        s = argv[i]

        if (s == "--user"):
            i += 1
            if (i > len(argv)):
                usage()
            g_user = argv[i]
        elif (s == "--pass"):
            i += 1
            if (i > len(argv)):
                usage()
            g_pass = argv[i]
        elif (s == "--sprint"):
            i += 1
            if (i > len(argv)):
                usage()
            g_sprint = argv[i]
        elif (s == "--csv"):
            g_csv = True
        elif (s == "--verbose"):
            g_verbose = True
        elif (s == "-h" or s == "--h" or s == "-help" or s == "--help"):
            usage()
        else:
            unknown_arg(s)

        i += 1

    if (g_user is None):
        print "ERROR: user is not specified"
        usage()

    if (g_pass is None):
        print "ERROR: pass is not specified"
        usage()

    if (g_sprint is None):
        print "ERROR: sprint is not specified"
        usage()


def main(argv):
    """
    Main program.

    @return: none
    """
    global g_script_name

    g_script_name = os.path.basename(argv[0])
    parse_config_file()
    parse_args(argv)

    url = 'https://0xdata.atlassian.net/rest/api/2/search?jql=sprint="' + urllib.quote(g_sprint) + '"&maxResults=1000'
    r = requests.get(url, auth=(g_user, g_pass))
    if (r.status_code != 200):
        print("ERROR: status code is " + str(r.status_code))
        sys.exit(1)
    j = r.json()
    issues = j[u'issues']
    pm = PeopleManager()
    for issue in issues:
        pm.add(issue)

    pm.emit()
    print("")

if __name__ == "__main__":
    main(sys.argv)
