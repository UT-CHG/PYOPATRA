// Method Definitions for List

template <class T>
void List<T>::push(ILLNode<T> &node) {
    if (tail) {
        node.insert_after(*tail);
        tail = &node;
    } else {
        head = &node;
        tail = &node;
    }
    length++;
}

template <class T>
void List<T>::remove(ILLNode<T> &node) {
    node.remove();
    length--;
}

//template <class T>
//ILLNode<T>* List<T>::pop_current() {
//    ILLNode<T>* temp = current;
//    current = temp->prev;
//    temp->remove();
//    return temp;
//}

//template <class T>
//void List<T>::advance() {
//    current = current->next;
//}