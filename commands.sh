
run(){
        python $PROJECT_PATH/scripts/$1.py "$@""
}

show(){
	python $PROJECT_PATH/scripts/visualisation.py $@"
}
