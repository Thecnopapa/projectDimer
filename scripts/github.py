import os
from utilities import *


def automatic_push_to_branch(target="auto", force = True, message="automated commit"):
  sprint("Automatic push to GitHub")
  print1("Target branch:", target)
  print1("Message:", message)
  import subprocess
  if target == "main":
    sprint("Are you crazy? don't push to main automatically!")
    return
  if len(target) == 0:
    sprint("Missing target")
    return

  #git_switch = ["git", "switch", "-C", target]
  git_switch = ["git", "checkout", "-b", target]
  git_add = ["git", "add", "*"]
  git_commit = ["git", "commit", "-a", "-m", message]
  git_push = ["git", "push", "-u", "origin/{} ", target, "--force"]
  try:
    print1("Switching branch")
    print2(" ".join(git_switch))
    subprocess.run(git_switch)
  except:
    print1("Switch failed")
    return
  
  try:
    print1("Adding changes")
    print2(" ".join(git_add))
    subprocess.run(git_add)
  except:
    print1("Addition failed")
    return

  try:
    print1("Committing changes")
    print2(" ".join(git_commit))
    subprocess.run(git_commit)
  except:
    print1("Commit failed")
    return

  try:
    print1("Pushing changes")
    print2(" ".join(git_push))
    subprocess.run(git_push)
  except:
    print1("Push failed")
    return

  print1("Automated git push successfull")


if __name__ == "__main__":
  automatic_push_to_branch(target="test")
  print("Done")