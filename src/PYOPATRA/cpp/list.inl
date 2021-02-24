// Method Definitions for List

template <class T>
void List<T>::push(ILLNode<T> &node) {
    node.insert_after(*tail);
    length++;
}

template <class T>
void List<T>::remove(ILLNode<T> &node) {
    current = node->prev;
    node.remove();
    length--;
}

ILLNode<T>* pop_current() {
    ILLNode<T>* temp = current;
    current = temp->prev;
    temp.remove();
    return temp;
}

template <class T>
void List<T>::advance() {
    current = current->next;
}