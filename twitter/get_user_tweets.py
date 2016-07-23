#!/usr/bin/env python

from config import *
import tweepy
import argparse

parser = argparse.ArgumentParser(description='Live Twitter streaming.')
parser.add_argument('-u','--user', default='SciPyTip', help='Get tweets from user (default: SciPyTip).', required=False)
parser.add_argument('-t','--topics', nargs='+', help='Topics that you want to follow.', required=False)
args = vars(parser.parse_args())

# Authentification handler and API
auth = tweepy.OAuthHandler(consumer_key, consumer_secret)
auth.set_access_token(access_token, access_secret)
api = tweepy.API(auth)


for status in tweepy.Cursor(api.user_timeline, screen_name=args['user']).items(200):
	for topic in args['topics']:
		if topic.lower() in status.text.lower():
			print("\n{}\n".format(status.text))
