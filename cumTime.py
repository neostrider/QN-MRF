file = open('timeTakenTRN.txt', 'r')

target = open('cumTimeTRN.txt', 'w')

timeTakenStr = file.readlines()

timeTaken = [float(x) for x in timeTakenStr]

cumTime = 0

index = 0

for x in timeTaken:
 cumTime = cumTime + x

 if ((index % 5) == 0):
  target.write(str(cumTime))
  target.write('\n')

 index = index + 1
