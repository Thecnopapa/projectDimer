
export PROJECT_PATH=$PROJECT_PATH

run(){
	python $PROJECT_PATH/scripts/$1.py "$@"
}

show(){
	python $PROJECT_PATH/scripts/visualisation.py "$@"
}

get_weights(){
	python $PROJECT_PATH/scripts/_get_weights.py "$@"
}