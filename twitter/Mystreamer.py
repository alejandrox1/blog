#!/usr/bin/env python

import argparse
from config import *
import tweepy
import json

### COMMAND LINE ARGUMENTS
parser = argparse.ArgumentParser(description='Live Twitter streaming.')
parser.add_argument('-o','--output', default='mytweets.json', help='Tweets will be written to this file (default: mytweets.json).', required=False)
parser.add_argument('-t','--topics', default='#Python', nargs='+', help='Trending topics that you want to follow (default: #Python).', required=False)
args = vars(parser.parse_args())
print(args)

### AUTHENTIFICATION HANDLER AND API
auth = tweepy.OAuthHandler(consumer_key, consumer_secret)
auth.set_access_token(access_token, access_secret)
api = tweepy.API(auth)

### PROCESS TWEETS
class MyListener(tweepy.streaming.StreamListener):
 
	def on_data(self, data):
		try:
			with open(args['output'], 'a') as f:
				f.write(data)
				# Parsing the tweet
				tweet = json.loads(data)
				print("{}".format(tweet['text']))       # If you want to see something else see the file 'TWEET_STRUCT'
				return True
		except BaseException as e:
			print("Error on_data: {}".format(e))
		return True

	def on_error(self, status):
		print(status)
		return True

	def save_bigram(self, bigram, most_frq):
		with open(args['screen'], 'a') as f:
			f.write('{}\n'.format(most_frq))
			for i in bigram:
				f.write('{}\n'.format(i))
			f.write('\n\n\n')
		return


if __name__=='__main__':
	twitter_stream = tweepy.Stream(auth, MyListener())
	twitter_stream.filter(track=args['topics'])

