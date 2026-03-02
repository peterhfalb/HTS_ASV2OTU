import sys
import traceback

def Die(Msg):
	print("",f=sys.stderr)
	print("",f=sys.stderr)

	traceback.print_stack()
	s = ""
	for i in range(0, len(sys.argv)):
		if i > 0:
			s += " "
		s += sys.argv[i]
	print(s, f=sys.stderr)
	print("**Error**", msg, f=sys.stderr)
	print("",f=sys.stderr)
	print("",f=sys.stderr)
	sys.exit(1)
	print("NOTHERE!!")

def Warning(Msg):
	print(Msg,  f=sys.stderr)
	print(sys.argv, f=sys.stderr)
	print("**Warning**",Msg, f=sys.stderr)
